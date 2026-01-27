import gzip
import os
import re
from logging import getLogger
from typing import Any, Dict, List, Optional

from ._parsers import _Parser

logger = getLogger(__name__)


class AnnotatedCnvVcfParser(_Parser):
    """
    Parser output is shaped to be directly mappable to the Django model fields.

    - Only returns FILTER == PASS records
    - VEP CSQ fields are treated as annotations (not caller_info)
    - caller_info contains raw caller-level fields (QUAL/ALT/raw INFO minus CSQ/raw sample)
    """

    _EXON_RE = re.compile(r"exon[s]?\s+(\d+)(?:-(\d+))?", re.IGNORECASE)

    def parse(self, vcf_path: str) -> List[Dict[str, Any]]:
        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        vcf_text = self._read_text_maybe_gzip(vcf_path)
        all_lines = [line for line in vcf_text.split("\n") if line]

        # Keep meta header lines so we can parse CSQ "Format:" definition
        meta_lines = [l for l in all_lines if l.startswith("##")]
        non_meta_lines = [l for l in all_lines if not l.startswith("##")]

        header_line = next((line for line in non_meta_lines if line.startswith("#CHROM")), None)
        if not header_line:
            raise ValueError("No header line found in VCF")

        headers = header_line.split("\t")
        idx = self._index_map(headers)

        csq_fields = self._parse_csq_format_fields(meta_lines)

        data_lines = [line for line in non_meta_lines if not line.startswith("#")]
        results: List[Dict[str, Any]] = []

        for line in data_lines:
            fields = line.split("\t")

            filter_status = self._safe_get(fields, idx["FILTER"], default="")
            if filter_status != "PASS":
                continue

            chrom = self._safe_get(fields, idx["#CHROM"], default="")
            pos = self._safe_int(self._safe_get(fields, idx["POS"]), default=0)

            vcf_id = self._safe_get(fields, idx["ID"], default="")
            alt = self._safe_get(fields, idx["ALT"], default="")
            qual = self._safe_float(
                self._safe_get(fields, idx["QUAL"]), default=float("nan")
            )

            info_string = self._safe_get(fields, idx["INFO"], default="")
            format_string = self._safe_get(fields, idx["FORMAT"], default="")
            sample_data = self._safe_get(fields, idx["SAMPLE"], default="")

            info = self._parse_info(info_string)
            sample_info = self._parse_sample(format_string, sample_data)

            # Caller-level SV fields
            sv_type_info = info.get("SVTYPE") or ""
            sv_len = self._safe_int(info.get("SVLEN"), default=0)
            end_pos = self._safe_int(info.get("END"), default=0)
            pos_range = [pos, end_pos]

            # ---- Parse VEP CSQ into annotation fields (flattened, no "vep_csq" key) ----
            vep_csq: Dict[str, Any] = {}
            csq_value = info.get("CSQ")
            if csq_value:
                first_csq = str(csq_value).split(",")[0]
                parts = first_csq.split("|")

                if csq_fields:
                    vep_csq = {
                        csq_fields[i]: (parts[i] if i < len(parts) else "")
                        for i in range(len(csq_fields))
                    }
                else:
                    # If header is missing, still keep the raw CSQ string as an annotation value
                    vep_csq = {"CSQ": first_csq}

            def _csq_get(*names: str) -> str:
                for n in names:
                    v = vep_csq.get(n)
                    if v not in (None, "", "."):
                        return str(v)
                return ""

            consequence = _csq_get("Consequence")
            impact = _csq_get("IMPACT")
            gene_symbol = _csq_get("SYMBOL")

            # Backward-compatible fallback only when we don't have CSQ header names
            if csq_value and (not consequence or not impact or not gene_symbol) and not csq_fields:
                parts = str(csq_value).split(",")[0].split("|")
                consequence = consequence or (parts[1] if len(parts) > 1 else "")
                impact = impact or (parts[2] if len(parts) > 2 else "")
                gene_symbol = gene_symbol or (parts[3] if len(parts) > 3 else "")

            clinical_significance = _csq_get("CLIN_SIG", "Clinical_significance", "clin_sig")
            cadd = self._safe_float(
                _csq_get("CADD", "CADD_PHRED", "CADD_PHRED_score", "CADD_score"),
                default=None,
            )
            p_haplo = self._safe_float(_csq_get("pHaplo", "p_haplo"), default=None)
            p_triplo = self._safe_float(_csq_get("pTriplo", "p_triplo"), default=None)

            # gnomAD SV AF is annotation (custom INFO)
            gnomad_sv_af: Optional[float] = None
            if info.get("gnomAD_SV_AF"):
                af_values = str(info["gnomAD_SV_AF"]).split("&")
                if af_values and af_values[0] not in ("", "."):
                    gnomad_sv_af = self._safe_float(af_values[0], default=None)

            # ---- sample-derived fields ----
            copy_number = self._safe_int(sample_info.get("CN"), default=None)

            quality_score = None
            if sample_info.get("SM") not in (None, "", "."):
                sm = self._safe_float(sample_info.get("SM"), default=None)
                quality_score = (sm * 100) if sm is not None else None

            # ---- derive CNV type (VEP-side interpretation) ----
            cnv_type = "Unknown"
            if sv_type_info == "CNV":
                if alt == "<DEL>":
                    cnv_type = "Deletion"
                elif alt == "<DUP>":
                    cnv_type = "Duplication"
            if cnv_type == "Duplication" and copy_number is not None and copy_number >= 5:
                cnv_type = "Amplification"

            # ---- exon labeling ----
            exons = ""
            start_exon = None
            end_exon = None
            total_exons = None

            m = self._EXON_RE.search(consequence or "")
            if m:
                start_exon = int(m.group(1))
                end_exon = int(m.group(2)) if m.group(2) else start_exon
                exons = (
                    f"Exon{'s' if end_exon != start_exon else ''} "
                    f"{start_exon}{('-' + str(end_exon)) if end_exon != start_exon else ''}"
                )
            elif "transcript_ablation" in (consequence or ""):
                exons = "Complete gene"
            elif "intron" in (consequence or ""):
                exons = "Intronic"

            # Flatten all VEP CSQ fields directly into annotations
            # and avoid duplicates like:
            #   SYMBOL vs gene_symbol, Consequence vs consequence, IMPACT vs impact,
            #   CLIN_SIG vs clinical_significance, CADD_PHRED vs cadd, gnomAD_SV_AF vs gnomad_sv_af
            annotations: Dict[str, Any] = dict(vep_csq)

            def _is_empty(v: Any) -> bool:
                return v in (None, "", ".")

            def _set_if_empty(key: str, value: Any) -> None:
                if key not in annotations or _is_empty(annotations.get(key)):
                    if not _is_empty(value):
                        annotations[key] = value

            # Ensure key VEP fields are populated (prefer CSQ-named fields)
            _set_if_empty("SYMBOL", gene_symbol or "Unknown")
            _set_if_empty("Consequence", consequence)
            _set_if_empty("IMPACT", impact)
            _set_if_empty("CLIN_SIG", clinical_significance or "Unknown")

            # Prefer storing CADD in the VEP-named field when available
            _set_if_empty("CADD_PHRED", cadd)

            # Prefer gnomAD SV AF under the VEP/custom-named field
            _set_if_empty("gnomAD_SV_AF", gnomad_sv_af)

            # Add derived (non-duplicate) annotation fields
            annotations["cnv_type"] = cnv_type
            annotations["exons"] = exons or "N/A"
            annotations["start_exon"] = start_exon
            annotations["end_exon"] = end_exon
            annotations["total_exons"] = total_exons

            # Human-friendly description derived from consequence
            consequence_for_desc = annotations.get("Consequence") or consequence or ""
            annotations["description"] = str(consequence_for_desc).replace("_", " ") if consequence_for_desc else "N/A"

            # Keep pHaplo/pTriplo as explicit annotation values (not typically present in CSQ headers)
            annotations["pHaplo"] = p_haplo
            annotations["pTriplo"] = p_triplo

            # caller_info should not include CSQ (VEP annotation)
            caller_info_info = dict(info)
            caller_info_info.pop("CSQ", None)

            caller_info: Dict[str, Any] = {
                "vcf_id": vcf_id,
                "alt": alt,
                "qual": qual,
                "svlen": sv_len,
                "end": end_pos,
                "quality_score": quality_score,
                "raw_info": caller_info_info,
                "raw_sample": sample_info,
            }

            # IMPORTANT: sv_type field is set from VEP-derived cnv_type (per request)
            results.append(
                {
                    "chrom": chrom,
                    "pos_range": pos_range,
                    "sv_type": cnv_type or None,
                    "copy_number": copy_number,
                    "filter_status": filter_status,
                    "annotations": annotations,
                    "caller_info": caller_info,
                }
            )

        return results

    def _read_text_maybe_gzip(self, path: str) -> str:
        if path.endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                return f.read()
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read()

    def _index_map(self, headers: List[str]) -> Dict[str, int]:
        def _idx(name: str) -> int:
            try:
                return headers.index(name)
            except ValueError:
                return -1

        format_idx = _idx("FORMAT")
        sample_idx = (format_idx + 1) if format_idx != -1 else -1

        return {
            "#CHROM": _idx("#CHROM"),
            "POS": _idx("POS"),
            "ID": _idx("ID"),
            "ALT": _idx("ALT"),
            "QUAL": _idx("QUAL"),
            "FILTER": _idx("FILTER"),
            "INFO": _idx("INFO"),
            "FORMAT": format_idx,
            "SAMPLE": sample_idx,
        }

    def _parse_info(self, info_string: str) -> Dict[str, Any]:
        info: Dict[str, Any] = {}
        for item in (info_string or "").split(";"):
            if not item:
                continue
            if "=" in item:
                k, v = item.split("=", 1)
                info[k] = v
            else:
                info[item] = None
        return info

    def _parse_sample(self, format_string: str, sample_data: str) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        if not format_string or not sample_data:
            return out
        keys = format_string.split(":")
        vals = sample_data.split(":")
        for i, k in enumerate(keys):
            out[k] = vals[i] if i < len(vals) else None
        return out

    def _safe_get(self, fields: List[str], index: int, default: Any = None) -> Any:
        if index is None or index < 0 or index >= len(fields):
            return default
        return fields[index]

    def _safe_int(self, value: Any, default: Any = 0) -> Any:
        if value in (None, "", "."):
            return default
        try:
            return int(value)
        except (TypeError, ValueError):
            return default

    def _safe_float(self, value: Any, default: Any = 0.0) -> Any:
        if value in (None, "", "."):
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    def _parse_csq_format_fields(self, meta_lines: List[str]) -> Optional[List[str]]:
        """
        Extract VEP CSQ field names from the VCF header line like:
          ##INFO=<ID=CSQ,...,Format: Allele|Consequence|IMPACT|SYMBOL|...>
        """
        csq_line = next((l for l in meta_lines if l.startswith("##INFO=<ID=CSQ")), None)
        if not csq_line or "Format: " not in csq_line:
            return None
        try:
            return csq_line.split("Format: ")[1].rstrip('">\n').split("|")
        except Exception as e:
            logger.warning(f"Failed to parse CSQ format header: {e}")
            return None