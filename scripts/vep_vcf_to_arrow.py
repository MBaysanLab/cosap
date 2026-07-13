#!/usr/bin/env python3
"""
Standalone script: read an Ensembl VEP-annotated VCF(.gz) and write it to an
Apache Arrow IPC Stream file using Polars for maximum speed.

Usage:
    python vep_vcf_to_arrow.py input.vcf.gz output.arrow
    python vep_vcf_to_arrow.py input.vcf.gz output.arrow --batch-size 500000
    python vep_vcf_to_arrow.py input.vcf.gz output.arrow --pick-only
    python vep_vcf_to_arrow.py input.vcf.gz output.arrow --pick-only --apply-transforms
    python vep_vcf_to_arrow.py input.vcf.gz output.arrow --compression zstd

Design goals:
    - Zero dependency on the cosap package.
    - Stream the file line-by-line so peak memory stays ≈ batch_size rows,
      not the full file.
    - All heavy lifting (CSQ explode, type casting, column rename) is done
      inside Polars expressions → benefits from its multithreaded engine.
    - The Arrow IPC Stream format is used (no footer with per-batch block
      offsets) so the writer holds O(1) memory regardless of batch count.
      Read back with ``pyarrow.ipc.open_stream`` or ``polars.read_ipc_stream``.

Requirements:
    pip install polars pyarrow tqdm
"""

from __future__ import annotations

import argparse
import gc
import gzip
import re
import sys
import time
from pathlib import Path
from typing import Optional

import polars as pl
import pyarrow as pa
import pyarrow.ipc as ipc
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Column mapping: internal name → VEP CSQ field name
# Mirrors cosap/workflows/variant_annotation/_config.py
# ---------------------------------------------------------------------------
VEP_COLUMN_MAPPING: dict[str, str] = {
    "variant_id": "Uploaded_variation",
    "chr": "#CHROM",
    "pos": "POS",
    "ref": "REF",
    "alt": "ALT",
    "gene_id": "Gene",
    "gene_symbol": "SYMBOL",
    "feature": "Feature",
    "feature_type": "Feature_type",
    "consequence": "Consequence",
    "cdna_position": "cDNA_position",
    "cds_position": "CDS_position",
    "rs_id": "Existing_variation",
    "impact": "IMPACT",
    "hgvsc": "HGVSc",
    "hgvsp": "HGVSp",
    "vep_pick": "PICK",
    "mane": "MANE",
    "canonical": "CANONICAL",
    "gnomadg_af": "gnomADg_AF",
    "gnomadg_afr_af": "gnomADg_AFR_AF",
    "gnomadg_ami_af": "gnomADg_AMI_AF",
    "gnomadg_amr_af": "gnomADg_AMR_AF",
    "gnomadg_asj_af": "gnomADg_ASJ_AF",
    "gnomadg_eas_af": "gnomADg_EAS_AF",
    "gnomadg_fin_af": "gnomADg_FIN_AF",
    "gnomadg_mid_af": "gnomADg_MID_AF",
    "gnomadg_nfe_af": "gnomADg_NFE_AF",
    "gnomadg_oth_af": "gnomADg_OTH_AF",
    "gnomadg_sas_af": "gnomADg_SAS_AF",
    "max_af": "MAX_AF",
    "aminoacid_change": "Amino_acids",
    "sift_score": "SIFT",
    "polyphen_score": "PolyPhen",
    "clinical_significance": "CLIN_SIG",
    "alphamissense_class": "am_class",
    "alphamissense_pathogenicity": "am_pathogenicity",
    "revel_score": "REVEL",
    "spliceai_pred_dp_ag": "SpliceAI_pred_DP_AG",
    "spliceai_pred_dp_al": "SpliceAI_pred_DP_AL",
    "spliceai_pred_dp_dg": "SpliceAI_pred_DP_DG",
    "spliceai_pred_dp_dl": "SpliceAI_pred_DP_DL",
    "spliceai_pred_ds_ag": "SpliceAI_pred_DS_AG",
    "spliceai_pred_ds_al": "SpliceAI_pred_DS_AL",
    "spliceai_pred_ds_dg": "SpliceAI_pred_DS_DG",
    "spliceai_pred_ds_dl": "SpliceAI_pred_DS_DL",
    "spliceai_pred_symbol": "SpliceAI_pred_SYMBOL",
    "cadd_phred": "CADD_PHRED",
    "cadd_raw": "CADD_RAW",
    "bayesdel_addaf": "BayesDel",
    "primateai": "PrimateAI",
    "pubmed": "PUBMED",
    "intervar_classification": "InterVar_classification",
    "intervar_pvs1": "InterVar_PVS1",
    "intervar_ps": "InterVar_PS",
    "intervar_pm": "InterVar_PM",
    "intervar_pp": "InterVar_PP",
    "intervar_ba1": "InterVar_BA1",
    "intervar_bs": "InterVar_BS",
    "intervar_bp": "InterVar_BP",
}

# Reverse: VEP field name → our internal name
_VEP_TO_INTERNAL = {v: k for k, v in VEP_COLUMN_MAPPING.items()}

# Numeric columns that should be cast to Float64 when possible
FLOAT_COLUMNS = {
    "gnomadg_af",
    "gnomadg_afr_af",
    "gnomadg_ami_af",
    "gnomadg_amr_af",
    "gnomadg_asj_af",
    "gnomadg_eas_af",
    "gnomadg_fin_af",
    "gnomadg_mid_af",
    "gnomadg_nfe_af",
    "gnomadg_oth_af",
    "gnomadg_sas_af",
    "max_af",
    "revel_score",
    "alphamissense_pathogenicity",
    "spliceai_pred_dp_ag",
    "spliceai_pred_dp_al",
    "spliceai_pred_dp_dg",
    "spliceai_pred_dp_dl",
    "spliceai_pred_ds_ag",
    "spliceai_pred_ds_al",
    "spliceai_pred_ds_dg",
    "spliceai_pred_ds_dl",
    "cadd_phred",
    "cadd_raw",
    "bayesdel_addaf",
    "primateai",
}

INT_COLUMNS = {
    "pos",
}

SPLICEAI_DS_COLUMNS = [
    "spliceai_pred_ds_ag",
    "spliceai_pred_ds_al",
    "spliceai_pred_ds_dg",
    "spliceai_pred_ds_dl",
]

SPLICEAI_THRESHOLDS = {"high_impact": 0.8, "moderate_impact": 0.5, "low_impact": 0.2}
REVEL_THRESHOLDS = {"pathogenic": 0.75, "likely_pathogenic": 0.5, "uncertain": 0.25}

# Strings that should be treated as null
NULL_STRINGS = {"", ".", "-", "None", "none", "NA", "N/A", "na"}

# Valid compression codecs for Arrow IPC
VALID_COMPRESSIONS = {"zstd", "lz4", "uncompressed"}


# ---------------------------------------------------------------------------
# VCF header parsing
# ---------------------------------------------------------------------------
def _open_vcf(path: str):
    """Return a text-mode file handle (transparently decompresses .gz)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _parse_csq_fields(meta_lines: list[str]) -> list[str]:
    """
    Extract CSQ sub-field names from the VCF ``##INFO=<ID=CSQ,...>`` header.

    Returns e.g. ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', ...].
    """
    for line in meta_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            if "Format: " not in line:
                raise ValueError(
                    "Found CSQ INFO header but it lacks a 'Format: ' definition"
                )
            return line.split("Format: ")[1].rstrip('">\n').split("|")
    raise ValueError(
        "No ##INFO=<ID=CSQ,...> header found — is this a VEP-annotated VCF?"
    )


def _parse_info_ids(meta_lines: list[str]) -> list[str]:
    """Extract INFO IDs from ``##INFO=<ID=...>`` header lines."""
    info_ids: list[str] = []
    for line in meta_lines:
        if not line.startswith("##INFO=<ID="):
            continue
        info_id = line[len("##INFO=<ID=") :].split(",", 1)[0].strip()
        # CSQ is parsed separately and must not be retained as a generic INFO
        # column, otherwise the large raw CSQ payload gets duplicated.
        if info_id and info_id != "CSQ":
            info_ids.append(info_id)
    return info_ids


def _normalize_info_name(info_id: str) -> str:
    """Map a VCF INFO ID to a prefixed, collision-resistant column name."""
    normalized = re.sub(r"[^0-9A-Za-z]+", "_", info_id.strip()).strip("_")
    if not normalized:
        raise ValueError(f"Unable to normalize INFO field name: {info_id!r}")
    return f"info_{normalized.lower()}"


def _build_info_column_map(info_ids: list[str]) -> dict[str, str]:
    """Build a deterministic INFO ID -> output column mapping."""
    info_columns: dict[str, str] = {}
    seen: dict[str, str] = {}
    for info_id in info_ids:
        column_name = _normalize_info_name(info_id)
        existing = seen.get(column_name)
        if existing is not None and existing != info_id:
            raise ValueError(
                "INFO field name collision after normalization: "
                f"{existing!r} and {info_id!r} both map to {column_name!r}"
            )
        seen[column_name] = info_id
        info_columns[info_id] = column_name
    return info_columns


def _parse_vcf_header(path: str) -> tuple[list[str], list[str], list[str], dict[str, str]]:
    """
    Read meta-lines and the #CHROM header line.

    Returns (meta_lines, vcf_columns, csq_fields, info_columns).
    """
    meta_lines: list[str] = []
    vcf_columns: list[str] = []
    with _open_vcf(path) as fh:
        for line in fh:
            if line.startswith("##"):
                meta_lines.append(line.rstrip("\n"))
            elif line.startswith("#CHROM"):
                vcf_columns = line.rstrip("\n").split("\t")
                break
            else:
                break
    if not vcf_columns:
        raise ValueError("No #CHROM header line found in VCF")
    csq_fields = _parse_csq_fields(meta_lines)
    info_columns = _build_info_column_map(_parse_info_ids(meta_lines))
    return meta_lines, vcf_columns, csq_fields, info_columns


def _parse_info_values(info: str, info_columns: dict[str, str]) -> dict[str, Optional[str]]:
    """Parse all INFO values from a raw INFO string."""
    values = {column_name: None for column_name in info_columns.values()}
    if not info:
        return values

    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, "1"
        column_name = info_columns.get(key)
        if column_name is not None:
            values[column_name] = value
    return values


# ---------------------------------------------------------------------------
# Batch reader — yields lists of raw data-line strings
# ---------------------------------------------------------------------------
def _iter_data_lines(path: str, batch_size: int):
    """
    Yield batches of data lines (strings) from the VCF, skipping all header
    lines.  Each yielded list has at most *batch_size* elements.
    """
    buf: list[str] = []
    with _open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            buf.append(line.rstrip("\n"))
            if len(buf) >= batch_size:
                yield buf
                buf = []
    if buf:
        yield buf


# ---------------------------------------------------------------------------
# Core: parse a batch of lines into a Polars DataFrame
# ---------------------------------------------------------------------------
def _parse_batch(
    lines: list[str],
    vcf_columns: list[str],
    csq_fields: list[str],
    info_columns: dict[str, str],
    apply_transforms: bool,
    pick_only: bool,
) -> pl.DataFrame:
    """
    Convert raw VCF data lines into a fully-processed Polars DataFrame.

    Memory-conscious version: only materialises CSQ sub-fields that are
    actually referenced by ``VEP_COLUMN_MAPPING``. Unmapped CSQ fields are
    discarded during the split, never allocated as Polars columns. This is
    the single biggest peak-memory saving for VCFs with rich CSQ payloads
    (50+ fields × millions of exploded rows).

    Steps:
      1. Split each line on TAB → fixed VCF columns.
      2. Extract the CSQ value from the INFO column.
      3. Explode multi-transcript CSQ entries (comma-separated).
      4. Split each CSQ entry on '|', keeping only mapped sub-fields.
      5. Build variant_id.
      6. Rename VEP field names → internal names.
      7. Optionally apply transforms (type casts, score interpretations, etc.).
    """
    # --- Step 1: split into VCF columns -----------------------------------
    # We only need CHROM, POS, REF, ALT, and INFO (index 7 in standard VCF).
    chroms: list[str] = []
    positions: list[str] = []
    refs: list[str] = []
    alts: list[str] = []
    csq_raw_list: list[Optional[str]] = []
    info_value_lists = {column_name: [] for column_name in info_columns.values()}

    for line in lines:
        parts = line.split("\t", 8)  # split at most into 9 parts (fast)
        chroms.append(parts[0])
        positions.append(parts[1])
        refs.append(parts[3])
        alts.append(parts[4])
        # Extract CSQ from INFO field (index 7)
        info = parts[7] if len(parts) > 7 else ""
        csq_val = _extract_csq(info)
        csq_raw_list.append(csq_val)
        info_values = _parse_info_values(info, info_columns)
        for column_name, values in info_value_lists.items():
            values.append(info_values[column_name])

    # Build a base DataFrame with the VCF-level columns
    base_columns = {
        "#CHROM": chroms,
        "POS": positions,
        "REF": refs,
        "ALT": alts,
        "_csq_raw": csq_raw_list,
        **info_value_lists,
    }
    base = pl.DataFrame(base_columns)
    # Free the Python-side lists immediately — they hold ~half a GB of
    # Python string overhead per 200k-row batch.
    del chroms, positions, refs, alts, csq_raw_list, info_value_lists, base_columns

    # --- Step 2-4: explode CSQ -------------------------------------------
    # Determine which CSQ sub-fields we actually want to keep.
    # Anything not in VEP_COLUMN_MAPPING is dropped at split-time so it
    # never becomes a column. This is the key memory optimisation.
    wanted_vep_fields = set(VEP_COLUMN_MAPPING.values())
    keep_csq_indices = [
        (i, name) for i, name in enumerate(csq_fields) if name in wanted_vep_fields
    ]

    n_csq = len(csq_fields)

    # Explode comma-separated transcripts
    df = (
        base.lazy()
        .filter(pl.col("_csq_raw").is_not_null())
        .with_columns(pl.col("_csq_raw").str.split(",").alias("_csq_list"))
        .drop("_csq_raw")
        .explode("_csq_list")
    )
    del base

    # Split each CSQ entry into the individual sub-fields, but only
    # project out the ones we actually need.
    csq_split_expr = pl.col("_csq_list").str.split_exact("|", n_csq - 1)
    df = df.with_columns(
        [
            csq_split_expr.struct.field(f"field_{idx}").alias(name)
            for idx, name in keep_csq_indices
        ]
    ).drop("_csq_list")

    # Materialise. Use streaming so Polars can pipeline the explode +
    # projection without holding the full intermediate.
    try:
        df = df.collect(engine="streaming")
    except TypeError:
        # Older Polars versions
        df = df.collect(streaming=True)

    if df.is_empty():
        return df

    # --- Step 5: variant_id -----------------------------------------------
    df = df.with_columns(
        (
            pl.col("#CHROM")
            + pl.lit("_")
            + pl.col("POS")
            + pl.lit("_")
            + pl.col("REF")
            + pl.lit("_")
            + pl.col("ALT")
        ).alias("variant_id")
    )

    # Ensure variant_id has chr prefix and replace '/' → '_'
    df = df.with_columns(
        pl.when(pl.col("variant_id").str.starts_with("chr"))
        .then(pl.col("variant_id"))
        .otherwise(pl.lit("chr") + pl.col("variant_id"))
        .str.replace_all("/", "_")
        .alias("variant_id"),
    )

    # --- Step 6: rename VEP → internal ------------------------------------
    cols_present = set(df.columns)
    rename = {
        vep_name: internal_name
        for internal_name, vep_name in VEP_COLUMN_MAPPING.items()
        if vep_name in cols_present and vep_name != internal_name
    }
    if rename:
        df = df.rename(rename)

    # Keep only mapped columns (+ variant_id which is already named correctly)
    keep = [c for c in VEP_COLUMN_MAPPING.keys() if c in df.columns]
    info_keep = [c for c in info_columns.values() if c in df.columns]
    df = df.select(keep + info_keep)

    # --- Step 7: optional transforms (mirrors cosap post-processing) ------
    if apply_transforms:
        df = _apply_transforms(df, pick_only)
    elif pick_only and "vep_pick" in df.columns:
        df = df.filter(pl.col("vep_pick") == "1")

    return df


def _extract_csq(info: str) -> Optional[str]:
    """
    Fast extraction of the CSQ=... value from a VCF INFO string.

    INFO fields are semicolon-separated key=value pairs.  We want
    the value after ``CSQ=`` up to the next ``;`` (or end of string).
    """
    # Fast path: look for 'CSQ=' directly (avoids splitting the whole string)
    idx = info.find("CSQ=")
    if idx == -1:
        return None
    start = idx + 4  # skip 'CSQ='
    end = info.find(";", start)
    if end == -1:
        return info[start:]
    return info[start:end]


# ---------------------------------------------------------------------------
# Transforms — mirrors annotation_output_parsers.py & _config.py
# ---------------------------------------------------------------------------
def _apply_transforms(df: pl.DataFrame, pick_only: bool) -> pl.DataFrame:
    """Apply the same transforms the cosap VEPResultsProcessor applies."""
    cols = set(df.columns)

    # ---- null normalisation -----------------------------------------------
    # Replace known null sentinel strings with actual nulls
    str_cols = [
        c for c, dt in zip(df.columns, df.dtypes) if dt == pl.Utf8 or dt == pl.String
    ]
    if str_cols:
        df = df.with_columns(
            [
                pl.when(pl.col(c).is_in(list(NULL_STRINGS)))
                .then(pl.lit(None))
                .otherwise(pl.col(c))
                .alias(c)
                for c in str_cols
            ]
        )

    # ---- numeric casts ----------------------------------------------------
    for c in FLOAT_COLUMNS & cols:
        df = df.with_columns(pl.col(c).cast(pl.Float64, strict=False))

    for c in INT_COLUMNS & cols:
        df = df.with_columns(pl.col(c).cast(pl.Int64, strict=False))

    # ---- boolean columns --------------------------------------------------
    if "canonical" in cols:
        df = df.with_columns(
            (pl.col("canonical").str.to_lowercase() == "yes").alias("canonical")
        )

    if "vep_pick" in cols:
        df = df.with_columns((pl.col("vep_pick") == "1").alias("vep_pick"))

    # ---- filter to PICK == true if requested ------------------------------
    if pick_only and "vep_pick" in cols:
        df = df.filter(pl.col("vep_pick"))

    # ---- consequence: take most severe (first before '&') -----------------
    if "consequence" in cols:
        df = df.with_columns(
            pl.col("consequence").str.split("&").list.first().alias("consequence")
        )

    # ---- SpliceAI max score + interpretation ------------------------------
    avail_splice = [c for c in SPLICEAI_DS_COLUMNS if c in cols]
    if avail_splice:
        df = df.with_columns(
            pl.max_horizontal(*[pl.col(c) for c in avail_splice]).alias(
                "spliceai_max_score"
            )
        )
        t = SPLICEAI_THRESHOLDS
        df = df.with_columns(
            pl.when(pl.col("spliceai_max_score") >= t["high_impact"])
            .then(pl.lit("high_impact"))
            .when(pl.col("spliceai_max_score") >= t["moderate_impact"])
            .then(pl.lit("moderate_impact"))
            .when(pl.col("spliceai_max_score") >= t["low_impact"])
            .then(pl.lit("low_impact"))
            .otherwise(pl.lit("minimal_impact"))
            .alias("spliceai_interpretation")
        )

    # ---- REVEL interpretation ---------------------------------------------
    if "revel_score" in cols:
        t = REVEL_THRESHOLDS
        df = df.with_columns(
            pl.when(pl.col("revel_score") >= t["pathogenic"])
            .then(pl.lit("pathogenic"))
            .when(pl.col("revel_score") >= t["likely_pathogenic"])
            .then(pl.lit("likely_pathogenic"))
            .when(pl.col("revel_score") >= t["uncertain"])
            .then(pl.lit("uncertain"))
            .otherwise(pl.lit("benign"))
            .alias("revel_interpretation")
        )

    # ---- clinical_significance → list of strings --------------------------
    if "clinical_significance" in cols:
        df = df.with_columns(
            pl.col("clinical_significance")
            .str.strip_chars()
            .str.split(",")
            .alias("clinical_significance")
        )

    return df


# ---------------------------------------------------------------------------
# Schema alignment — make sure every batch has the same Arrow schema
# ---------------------------------------------------------------------------
def _align_to_schema(table: pa.Table, schema: pa.Schema) -> pa.Table:
    """
    Re-order / pad *table* so it exactly matches *schema*.

    - Columns in *schema* but missing in *table* get filled with nulls.
    - Extra columns in *table* that are not in *schema* are dropped.
    - Column types are cast when safe.
    """
    arrays = []
    for field in schema:
        if field.name in table.column_names:
            col = table.column(field.name)
            if col.type != field.type:
                try:
                    col = col.cast(field.type, safe=False)
                except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                    pass  # keep as-is, Arrow IPC handles this gracefully
            arrays.append(col)
        else:
            arrays.append(pa.nulls(len(table), type=field.type))
    return pa.table(arrays, schema=schema)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def convert(
    input_path: str,
    output_path: str,
    batch_size: int = 50_000,
    pick_only: bool = False,
    apply_transforms: bool = False,
    compression: str = "uncompressed",
) -> None:
    """
    Read *input_path* (VEP VCF / VCF.gz) and write *output_path* (Arrow IPC
    Stream).

    Parameters
    ----------
    input_path : str
        Path to the VEP-annotated VCF (plain or gzipped).
    output_path : str
        Path to the output Arrow IPC (.arrow) file.
    batch_size : int
        Number of VCF data lines to process per batch.  Higher values use more
        memory but benefit from Polars vectorisation.  200k is a good default.
    pick_only : bool
        If True, only keep rows where VEP PICK == 1 (one row per variant).
    apply_transforms : bool
        If True, apply the full cosap-style post-processing (type casts,
        score interpretations, null normalisation, column renaming).
        If False, all CSQ fields are kept as raw strings (fastest).
    compression : str
        Compression codec for Arrow IPC record batches.
        One of "zstd", "lz4", or "uncompressed" (default).
    """
    compression = compression.lower()
    if compression not in VALID_COMPRESSIONS:
        raise ValueError(
            f"Unsupported compression: {compression!r}. "
            f"Choose from: {', '.join(sorted(VALID_COMPRESSIONS))}"
        )

    t0 = time.perf_counter()

    # ---- parse header ----
    _meta, vcf_columns, csq_fields, info_columns = _parse_vcf_header(input_path)
    print(f"VCF columns : {vcf_columns[:8]}{'...' if len(vcf_columns) > 8 else ''}")
    print(f"CSQ fields  : {len(csq_fields)} ({', '.join(csq_fields[:6])}, ...)")
    print(
        "INFO fields : "
        f"{len(info_columns)} ({', '.join(list(info_columns.values())[:6])}"
        f"{', ...' if len(info_columns) > 6 else ''})"
    )
    if compression != "uncompressed":
        print(f"Compression : {compression}")

    # ---- IPC write options ----
    ipc_options = (
        ipc.IpcWriteOptions(compression=compression)
        if compression != "uncompressed"
        else None
    )

    # ---- streaming write ----
    writer: Optional[ipc.RecordBatchStreamWriter] = None
    schema: Optional[pa.Schema] = None
    total_rows = 0
    with tqdm(desc="Processing VCF batches", unit="batch") as pbar:
        for line_batch in _iter_data_lines(input_path, batch_size):
            df = _parse_batch(
                line_batch,
                vcf_columns,
                csq_fields,
                info_columns,
                apply_transforms,
                pick_only,
            )

            batch_rows = len(df)
            if df.is_empty():
                del df
                pbar.update(1)
                pbar.set_postfix_str(f"rows=0 skipped total_rows={total_rows:,}")
                continue

            arrow_table = df.to_arrow()
            del df

            if writer is None:
                schema = arrow_table.schema
                writer = ipc.new_stream(output_path, schema, options=ipc_options)

            # Align this batch to the first batch's schema (handles column
            # order differences and missing columns across batches).
            arrow_table = _align_to_schema(arrow_table, schema)
            writer.write_table(arrow_table)

            # Release per-batch buffers back to the OS to prevent the slow
            # RSS creep that PyArrow's memory pool otherwise causes when
            # processing hundreds of batches.
            del arrow_table
            pa.default_memory_pool().release_unused()
            gc.collect()

            total_rows += batch_rows
            pbar.update(1)
            pbar.set_postfix_str(
                f"rows={batch_rows:,} total_rows={total_rows:,}"
            )

    if writer is not None:
        writer.close()

    elapsed = time.perf_counter() - t0
    size_mb = (
        Path(output_path).stat().st_size / (1024 * 1024)
        if Path(output_path).exists()
        else 0
    )
    print(
        f"\nDone — {total_rows:,} rows written to {output_path} "
        f"({size_mb:.1f} MB, {compression}) in {elapsed:.2f}s"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Convert an Ensembl VEP-annotated VCF(.gz) to Apache Arrow IPC.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fast raw dump (no transforms, all transcripts):
  python vep_vcf_to_arrow.py sample.vcf.gz sample.arrow

  # Keep only PICK==1 rows, apply cosap-style transforms:
  python vep_vcf_to_arrow.py sample.vcf.gz sample.arrow --pick-only --apply-transforms

  # Zstandard-compressed output (typically 3-5x smaller):
  python vep_vcf_to_arrow.py sample.vcf.gz sample.arrow --compression zstd

  # Large file with custom batch size and LZ4 compression:
  python vep_vcf_to_arrow.py huge.vcf.gz huge.arrow --batch-size 500000 --compression lz4
""",
    )
    parser.add_argument("input", help="Input VCF or VCF.gz file (VEP-annotated)")
    parser.add_argument("output", help="Output Arrow IPC file path (.arrow)")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50_000,
        help=(
            "Number of VCF lines per processing batch (default: 50000). "
            "Each line can explode into many rows (one per VEP transcript), "
            "so the effective row count per batch is often 10-50x this value. "
            "Lower if you hit OOM; raise if you have plenty of RAM."
        ),
    )
    parser.add_argument(
        "--pick-only",
        action="store_true",
        help="Only keep rows where VEP PICK == 1 (one row per variant)",
    )
    parser.add_argument(
        "--apply-transforms",
        action="store_true",
        help=(
            "Apply cosap-style post-processing: null normalisation, numeric "
            "casts, score interpretations, clinical_significance as list"
        ),
    )
    parser.add_argument(
        "--compression",
        choices=sorted(VALID_COMPRESSIONS),
        default="uncompressed",
        help=(
            "Compression codec for Arrow IPC record batches. "
            "'zstd' gives best size/speed trade-off, 'lz4' is fastest "
            "decompression (default: uncompressed)"
        ),
    )
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    convert(
        input_path=args.input,
        output_path=args.output,
        batch_size=args.batch_size,
        pick_only=args.pick_only,
        apply_transforms=args.apply_transforms,
        compression=args.compression,
    )


if __name__ == "__main__":
    main()
