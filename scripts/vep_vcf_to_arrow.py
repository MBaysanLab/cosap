#!/usr/bin/env python3
"""
Standalone script: read an Ensembl VEP-annotated VCF(.gz) and write it to an
Apache Arrow IPC file using Polars for maximum speed.

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
    - The Arrow IPC file is written incrementally (one RecordBatch per batch)
      so the output is valid even for huge files.

Requirements:
    pip install polars pyarrow
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path
from typing import Optional

import polars as pl
import pyarrow as pa
import pyarrow.ipc as ipc

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


def _parse_vcf_header(path: str) -> tuple[list[str], list[str], list[str]]:
    """
    Read meta-lines and the #CHROM header line.

    Returns (meta_lines, vcf_columns, csq_fields).
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
    return meta_lines, vcf_columns, csq_fields


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
    apply_transforms: bool,
    pick_only: bool,
) -> pl.DataFrame:
    """
    Convert raw VCF data lines into a fully-processed Polars DataFrame.

    Steps:
      1. Split each line on TAB → fixed VCF columns.
      2. Extract the CSQ value from the INFO column.
      3. Explode multi-transcript CSQ entries (comma-separated).
      4. Split each CSQ entry on '|' into individual fields.
      5. Build variant_id.
      6. Rename VEP field names → internal names.
      7. Optionally apply transforms (type casts, score interpretations, etc.).
    """
    n_vcf_cols = len(vcf_columns)

    # --- Step 1: split into VCF columns -----------------------------------
    # We only need CHROM, POS, REF, ALT, and INFO (index 7 in standard VCF).
    # Doing a full split is cheap relative to CSQ parsing.
    chroms: list[str] = []
    positions: list[str] = []
    refs: list[str] = []
    alts: list[str] = []
    csq_raw_list: list[Optional[str]] = []

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

    # Build a base DataFrame with the VCF-level columns
    base = pl.DataFrame(
        {
            "#CHROM": chroms,
            "POS": positions,
            "REF": refs,
            "ALT": alts,
            "_csq_raw": csq_raw_list,
        }
    )

    # --- Step 2-4: explode CSQ -------------------------------------------
    # Explode comma-separated transcripts, then split on '|'
    df = (
        base.filter(pl.col("_csq_raw").is_not_null())
        .with_columns(pl.col("_csq_raw").str.split(",").alias("_csq_list"))
        .explode("_csq_list")
        .drop("_csq_raw")
    )

    if df.is_empty():
        return df.drop("_csq_list")

    # Split each CSQ entry into the individual sub-fields
    n_csq = len(csq_fields)
    csq_struct = df.select(pl.col("_csq_list").str.split_exact("|", n_csq - 1)).unnest(
        "_csq_list"
    )

    # Rename struct fields (field_0 .. field_N) → actual CSQ names
    rename_map = {f"field_{i}": csq_fields[i] for i in range(n_csq)}
    csq_struct = csq_struct.rename(rename_map)

    df = pl.concat([df.drop("_csq_list"), csq_struct], how="horizontal")

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

    # Ensure variant_id has chr prefix
    df = df.with_columns(
        pl.when(pl.col("variant_id").str.starts_with("chr"))
        .then(pl.col("variant_id"))
        .otherwise(pl.lit("chr") + pl.col("variant_id"))
        .alias("variant_id"),
    )

    # Also replace '/' → '_' in variant_id (matches cosap behaviour
    # for Uploaded_variation from VEP which uses chr_pos_ref/alt)
    df = df.with_columns(
        pl.col("variant_id").str.replace_all("/", "_").alias("variant_id")
    )

    # --- Step 6: rename VEP → internal ------------------------------------
    # Only rename columns that actually exist in the DataFrame
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
    # Also keep any CSQ columns that weren't in the mapping (user might want
    # to inspect the full VEP output).
    unmapped_csq = [
        c for c in df.columns if c not in keep and c not in {"_csq_raw", "_csq_list"}
    ]
    df = df.select(keep + unmapped_csq)

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
    batch_size: int = 200_000,
    pick_only: bool = False,
    apply_transforms: bool = False,
    compression: str = "uncompressed",
) -> None:
    """
    Read *input_path* (VEP VCF / VCF.gz) and write *output_path* (Arrow IPC).

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
    _meta, vcf_columns, csq_fields = _parse_vcf_header(input_path)
    print(f"VCF columns : {vcf_columns[:8]}{'...' if len(vcf_columns) > 8 else ''}")
    print(f"CSQ fields  : {len(csq_fields)} ({', '.join(csq_fields[:6])}, ...)")
    if compression != "uncompressed":
        print(f"Compression : {compression}")

    # ---- IPC write options ----
    ipc_options = (
        ipc.IpcWriteOptions(compression=compression)
        if compression != "uncompressed"
        else None
    )

    # ---- streaming write ----
    writer: Optional[ipc.RecordBatchFileWriter] = None
    schema: Optional[pa.Schema] = None
    total_rows = 0
    batch_num = 0

    for line_batch in _iter_data_lines(input_path, batch_size):
        batch_num += 1
        df = _parse_batch(
            line_batch, vcf_columns, csq_fields, apply_transforms, pick_only
        )

        if df.is_empty():
            print(f"  batch {batch_num}: 0 rows (all filtered / no CSQ), skipped")
            continue

        arrow_table = df.to_arrow()

        if writer is None:
            schema = arrow_table.schema
            writer = ipc.new_file(output_path, schema, options=ipc_options)

        # Align this batch to the first batch's schema (handles column
        # order differences and missing columns across batches).
        arrow_table = _align_to_schema(arrow_table, schema)
        writer.write_table(arrow_table)

        total_rows += len(df)
        print(
            f"  batch {batch_num}: {len(df):>9,} rows  (cumulative {total_rows:>12,})"
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
        default=200_000,
        help="Number of VCF lines per processing batch (default: 200000)",
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
