import gzip
import os
from logging import getLogger
from subprocess import run
from typing import Optional

import orjson
import pandas as pd
import pysam
from cyvcf2 import VCF
from pandas.api.types import is_numeric_dtype

logger = getLogger(__name__)


class VCFUtils:
    @staticmethod
    def read_vcf_into_df(path: str) -> pd.DataFrame:
        import io

        def file_open(path):
            if path.endswith(".gz"):
                return gzip.open(path, "rt")
            return open(path, "r")

        with file_open(path) as f:
            lines = [l for l in f if not l.startswith("##")]
            df = pd.read_csv(
                io.StringIO("".join(lines)),
                dtype={
                    "#CHROM": str,
                    "POS": int,
                    "ID": str,
                    "REF": str,
                    "ALT": str,
                    "QUAL": str,
                    "FILTER": str,
                    "INFO": str,
                },
                sep="\t",
                index_col=False,
            ).rename(columns={"#CHROM": "CHROM"})

        df.reset_index(inplace=True)
        df.rename(columns={"index": "id"}, inplace=True)
        return df

    @staticmethod
    def remove_same_ref_alt_alleles(path: str):
        new_path = path + ".removed_same_ref_alt.vcf"
        if path.endswith(".gz"):
            # Keep header lines, and keep variants where:
            # - REF != ALT (column 4 != column 5)
            # - ALT is not <NON_REF>
            awk_command = rf"zcat {path} | awk -F '\t' '/^#/ || ($4 != $5 && $5 !~ /<NON_REF>/)' > {new_path}"
        else:
            awk_command = rf"awk -F '\t' '/^#/ || ($4 != $5 && $5 !~ /<NON_REF>/)' {path} > {new_path}"

        run(["bash", "-c", awk_command], check=True)
        return new_path

    @staticmethod
    def convert_vcf_to_tsv(path: str, caller_type: str = "gatk") -> str:
        # Gatk VariantsToTable gives error when ref and alt alleles are same, so remove those lines

        path = VCFUtils.remove_same_ref_alt_alleles(path)
        output_filename = path.replace(".vcf", ".tsv")

        if os.path.exists(output_filename):
            return output_filename

        # Run GATK VariantsToTable to convert vcf to tsv
        command = [
            "gatk",
            "VariantsToTable",
            "-V",
            path,
            "-F",
            "CHROM",
            "-F",
            "POS",
            "-F",
            "REF",
            "-F",
            "ALT",
            "-F",
            "QUAL",
            "-F",
            "FILTER",
            "-F",
            "STATUS",
            "-F",
            "AF",
            "-GF",
            "AF",
            "-GF",
            "AD",
            "-GF",
            "DP",
            "-GF",
            "GT",
            "-SMA",
            "-O",
            output_filename,
        ]

        # If caller is strelka, to calculate AF we need AU, CU, GU, TU
        if caller_type.lower() == "strelka":
            command.extend(["-GF", "AU", "-GF", "CU", "-GF", "GU", "-GF", "TU"])

        if caller_type.lower() == "varscan":
            command.extend(["-GF", "FREQ"])

        if caller_type.lower() == "varnet":
            command.extend(["-GF", "AO"])

        if caller_type.lower() == "somaticsniper":
            command.extend(["-GF", "DP4"])

        _ = run(command, check=True, capture_output=True, text=True)

        return output_filename

    @staticmethod
    def calculate_strelka_af(row: pd.Series) -> float:
        """
        Calculates AF from AU, CU, GU, TU columns.
        From https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
        """
        ref = row["REF"]
        alt = row["ALT"]
        ref_u = f"TUMOR.{ref}U"
        alt_u = f"TUMOR.{alt}U"
        ref_counts = row[ref_u]
        alt_counts = row[alt_u]
        tier1_refcounts = ref_counts.split(",")[0]
        tier1_altcounts = alt_counts.split(",")[1]
        af = float(tier1_altcounts) / (float(tier1_altcounts) + float(tier1_refcounts))
        return af

    @staticmethod
    def calculate_strelka_ad(row: pd.Series) -> float:
        """
        Calculates AD as sum of AU, CU, GU, TU columns.
        """
        return sum(
            [
                int(i)
                for i in row[["TUMOR.AU", "TUMOR.CU", "TUMOR.GU", "TUMOR.TU"]]
                .values[0]
                .split(",")
            ]
        )

    @staticmethod
    def calculate_somaticsniper_af_ad(row: pd.Series) -> float:
        """
        Calculates AF from DP and DP4 columns.
        """
        dp = row["TUMOR.DP"]
        dp4 = row["TUMOR.DP4"]
        dp4 = dp4.split(",")
        ref_counts = sum([int(i) for i in dp4[:2]])
        alt_counts = sum([int(i) for i in dp4[2:]])
        af = alt_counts / (alt_counts + ref_counts)
        return af, alt_counts

    @staticmethod
    def get_column_mappings(vcf_df, sample_name, caller_type):
        caller_type = caller_type.lower()

        if caller_type == "strelka":
            for i, row in vcf_df.iterrows():
                af = VCFUtils.calculate_strelka_af(row)
                vcf_df.at[i, "strelka_AF"] = af
                ad = VCFUtils.calculate_strelka_ad(row)
                vcf_df.at[i, "strelka_AD"] = ad

            return {
                "AF": "strelka_AF",
                "AD": "strelka_AD",
            }
        elif caller_type in ["haplotypecaller", "deepvariant"]:
            return {
                "AD": f"{sample_name}.AD",
                "DP": f"{sample_name}.DP",
                "GT": f"{sample_name}.GT",
            }
        elif caller_type == "varscan":
            vcf_df["varscan_AF"] = (
                vcf_df[f"{sample_name}.FREQ"].str.replace("%", "").astype(float) / 100
            )

            return {
                "AF": "varscan_AF",
                "AD": f"{sample_name}.AD",
                "DP": f"{sample_name}.DP",
                "GT": f"{sample_name}.GT",
            }
        elif caller_type in ["mutect", "mutect2", "vardict"]:
            return {
                "AF": f"{sample_name}.AF",
                "AD": f"{sample_name}.AD",
                "DP": f"{sample_name}.DP",
                "GT": f"{sample_name}.GT",
            }
        elif caller_type == "somaticsniper":
            for i, row in vcf_df.iterrows():
                af, ad = VCFUtils.calculate_somaticsniper_af_ad(row)
                vcf_df.at[i, "somaticsniper_AF"] = af
                vcf_df.at[i, "somaticsniper_AD"] = ad

            return {
                "AF": "somaticsniper_AF",
                "AD": "somaticsniper_AD",
                "DP": f"{sample_name}.DP",
                "GT": f"{sample_name}.GT",
            }

        elif caller_type == "varnet":
            return {
                "AF": "SAMPLE.AF",
                "AD": "SAMPLE.AO",
                "DP": "SAMPLE.DP",
                "GT": "SAMPLE.GT",
            }
        else:
            raise ValueError(f"Unknown tool type: {caller_type}")

    @staticmethod
    def convert_vcf_to_dict(
        path: str,
        caller_type: str = "mutect2",
        sample_name: str = None,
        filter: str = "PASS",
    ) -> list:
        """
        Returns list of variants as json objects.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

        if sample_name is None:
            sample_names = VCFUtils.extract_sample_names_from_vcf(path)
            if sample_names is None:
                raise ValueError(
                    "Sample name not provided and could not be extracted from VCF file."
                )
            elif len(sample_names) == 1:
                sample_name = sample_names[0]
            elif len(sample_names) > 1 and "tumor" in [s.lower() for s in sample_names]:
                sample_name = "tumor"
            elif len(sample_names) > 1:
                raise ValueError(f"Multiple samples found in VCF file: {sample_names}")
            else:
                raise ValueError("Sample name not provided and could not be extracted.")

        if caller_type is None:
            caller_type = VCFUtils.extract_variantcaller_from_vcf(path)
            if caller_type is None:
                raise ValueError(
                    "Caller type not provided and could not be extracted from VCF file."
                )

        variants_table = VCFUtils.convert_vcf_to_tsv(path, caller_type)
        vcf_df = pd.read_csv(variants_table, sep="\t")

        vcf_df = vcf_df[vcf_df["FILTER"] == filter]

        vcf_df.columns = vcf_df.columns.str.upper()
        sample_name = sample_name.upper()

        column_mappings = VCFUtils.get_column_mappings(vcf_df, sample_name, caller_type)

        columns_to_keep = {
            "CHROM": "CHROM",
            "POS": "POS",
            "REF": "REF",
            "ALT": "ALT",
            "AF": "AF",
            "AD": "AD",
            "DP": "DP",
            "GT": "GT",
            "STATUS": "STATUS",
        }

        for key, value in column_mappings.items():
            columns_to_keep[key] = value

        try:
            final_df = vcf_df.filter(columns_to_keep.values())
        except KeyError as e:
            raise KeyError(f"Column not found in VCF file: {e}")

        # Reverse the name mappings
        final_df.rename(
            columns={v: k for k, v in columns_to_keep.items()}, inplace=True
        )
        print(final_df.head())
        # Check if AD column exists and process it if needed
        if "AD" in final_df.columns and not is_numeric_dtype(final_df["AD"]):
            final_df["AD"] = final_df["AD"].apply(
                lambda x: x.split(",")[-1] if isinstance(x, str) and "," in x else x
            )

        return final_df.to_dict("records")

    @staticmethod
    def convert_vcf_to_json(
        path: str,
        caller_type: str = None,
        sample_name: str = None,
        filter: str = "PASS",
    ) -> list:
        """
        Returns list of variants as json objects.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

        if sample_name is None:
            sample_names = VCFUtils.extract_sample_names_from_vcf(path)
            if sample_names is None:
                raise ValueError(
                    "Sample name not provided and could not be extracted from VCF file."
                )
            elif len(sample_names) == 1:
                sample_name = sample_names[0]
            elif len(sample_names) > 1 and "tumor" in [s.lower() for s in sample_names]:
                sample_name = "tumor"
            elif len(sample_names) > 1:
                raise ValueError(f"Multiple samples found in VCF file: {sample_names}")
            else:
                raise ValueError("Sample name not provided and could not be extracted.")

        if caller_type is None:
            caller_type = VCFUtils.extract_variantcaller_from_vcf(path)
            if caller_type is None:
                caller_type = "haplotypecaller"

        variants_table = VCFUtils.convert_vcf_to_tsv(path, caller_type)
        vcf_df = pd.read_csv(variants_table, sep="\t")

        vcf_df = vcf_df[vcf_df["FILTER"] == filter]

        vcf_df.columns = vcf_df.columns.str.upper()
        sample_name = sample_name.upper()

        column_mappings = VCFUtils.get_column_mappings(vcf_df, sample_name, caller_type)

        columns_to_keep = {
            "CHROM": "CHROM",
            "POS": "POS",
            "REF": "REF",
            "ALT": "ALT",
            "AF": "AF",
            "AD": "AD",
            "DP": "DP",
            "GT": "GT",
            "STATUS": "STATUS",
        }

        for key, value in column_mappings.items():
            columns_to_keep[key] = value

        try:
            final_df = vcf_df.filter(columns_to_keep.values())
        except KeyError as e:
            raise KeyError(f"Column not found in VCF file: {e}")

        # Reverse the name mappings
        final_df.rename(
            columns={v: k for k, v in columns_to_keep.items()}, inplace=True
        )

        if "AD" in final_df.columns and not is_numeric_dtype(final_df["AD"]):
            final_df["AD"] = final_df["AD"].apply(
                lambda x: x.split(",")[-1] if "," in x else x
            )

        return final_df.to_json(orient="records")

    @staticmethod
    def convert_vcf_to_json_fast(
        path: str,
        caller_type: str = None,
        sample_name: str = None,
        filter: str = "PASS",
    ) -> str:
        """
        Returns list of variants as json objects using cyvcf2 for fast parsing.
        If no sample information is available, returns only CHROM, POS, REF, ALT, and INFO fields.
        """

        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

        vcf = VCF(path)

        # Check if VCF has sample information
        sample_names = vcf.samples
        has_samples = bool(sample_names)
        sample_idx = None

        if has_samples:
            # Extract sample name if not provided
            if sample_name is None:
                if len(sample_names) == 1:
                    sample_name = sample_names[0]
                elif len(sample_names) > 1 and "tumor" in [
                    s.lower() for s in sample_names
                ]:
                    sample_name = next(s for s in sample_names if s.lower() == "tumor")
                elif len(sample_names) > 1:
                    # Use first sample if multiple exist and no tumor
                    sample_name = sample_names[0]

            # Get sample index
            if sample_name:
                try:
                    sample_idx = vcf.samples.index(sample_name)
                except ValueError:
                    # Sample not found, proceed without sample info
                    has_samples = False
                    sample_idx = None

        # Extract caller type if not provided and we have samples
        if caller_type is None and has_samples:
            caller_type = VCFUtils.extract_variantcaller_from_vcf(path)
            if caller_type is None:
                caller_type = "haplotypecaller"

        if caller_type:
            caller_type = caller_type.lower()

        # Parse variants
        variants = []
        for variant in vcf:
            # Apply filter
            if filter:
                variant_filter = variant.FILTER
                if variant_filter is None:
                    # PASS or no filter
                    pass
                elif filter == "PASS" and variant_filter:
                    # Skip non-PASS variants when filter is "PASS"
                    continue
                elif variant_filter != filter:
                    continue

            # Extract basic fields (always available)
            chrom = str(variant.CHROM)
            pos = int(variant.POS)
            qual = float(variant.QUAL) if variant.QUAL is not None else None
            ref = str(variant.REF)
            alt = ",".join(variant.ALT) if variant.ALT else "."

            # If no sample information, return basic variant info only
            if not has_samples or sample_idx is None:
                variant_data = {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                }

                # Try to extract some INFO fields if available
                qual = variant.QUAL
                if qual is not None:
                    variant_data["QUAL"] = float(qual)

                # Extract common INFO fields
                info_fields_to_extract = ["AF", "DP", "AC", "AN", "STATUS"]
                for field in info_fields_to_extract:
                    try:
                        value = variant.INFO.get(field)
                        if value is not None:
                            if isinstance(value, tuple):
                                variant_data[field] = (
                                    value[0] if len(value) == 1 else list(value)
                                )
                            else:
                                variant_data[field] = value
                    except (KeyError, TypeError):
                        pass

                variants.append(variant_data)
                continue

            # Extract sample-specific fields
            status = (
                str(variant.INFO.get("STATUS", "."))
                if variant.INFO.get("STATUS") is not None
                else "."
            )

            gt = variant.gt_types[sample_idx] if hasattr(variant, "gt_types") else None

            dp_val = (
                variant.format("DP")[sample_idx][0] if "DP" in variant.FORMAT else None
            )
            dp = int(dp_val) if dp_val is not None else None

            # Calculate AF and AD based on caller type
            af = None
            ad = None

            if caller_type == "strelka":
                af, ad = VCFUtils._extract_strelka_af_ad_cyvcf2(
                    variant, sample_idx, ref, alt
                )
            elif caller_type == "somaticsniper":
                af, ad = VCFUtils._extract_somaticsniper_af_ad_cyvcf2(
                    variant, sample_idx
                )
            elif caller_type == "varscan":
                freq = (
                    variant.format("FREQ")[sample_idx][0]
                    if "FREQ" in variant.FORMAT
                    else None
                )
                if freq:
                    af = float(str(freq).replace("%", "")) / 100.0
                ad_values = (
                    variant.format("AD")[sample_idx] if "AD" in variant.FORMAT else None
                )
                ad = (
                    int(ad_values[-1])
                    if ad_values is not None and len(ad_values) > 1
                    else None
                )
            elif caller_type == "varnet":
                af_values = (
                    variant.format("AF")[sample_idx] if "AF" in variant.FORMAT else None
                )
                af = float(af_values[0]) if af_values is not None else None
                ad_values = (
                    variant.format("AO")[sample_idx] if "AO" in variant.FORMAT else None
                )
                ad = int(ad_values[0]) if ad_values is not None else None
            elif caller_type == "haplotypecaller":
                ad_values = (
                    variant.format("AD")[sample_idx] if "AD" in variant.FORMAT else None
                )
                ad = (
                    int(ad_values[-1])
                    if ad_values is not None and len(ad_values) > 1
                    else None
                )
                # Calculate AF from AD
                if ad_values is not None and len(ad_values) > 1:
                    total = sum(ad_values)
                    af = float(ad_values[-1]) / total if total > 0 else 0.0
            else:  # Default extraction for mutect, vardict, etc.
                af_values = (
                    variant.format("AF")[sample_idx] if "AF" in variant.FORMAT else None
                )
                af = float(af_values[0]) if af_values is not None else None
                ad_values = (
                    variant.format("AD")[sample_idx] if "AD" in variant.FORMAT else None
                )
                ad = (
                    int(ad_values[-1])
                    if ad_values is not None and len(ad_values) > 1
                    else None
                )

            # Format GT - explicitly convert to string
            gt_str = None
            if gt is not None:
                gt_bases = (
                    variant.gt_bases[sample_idx]
                    if hasattr(variant, "gt_bases")
                    else None
                )
                gt_str = str(gt_bases) if gt_bases else str(gt)

            # Ensure all numeric values are Python native types
            variants.append(
                {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "AF": float(af) if af is not None else None,
                    "AD": int(ad) if ad is not None else None,
                    "DP": dp,
                    "GT": gt_str,
                    "STATUS": status,
                }
            )

        vcf.close()

        # Use orjson to serialize
        return orjson.dumps(variants).decode("utf-8")

    @staticmethod
    def _extract_strelka_af_ad_cyvcf2(variant, sample_idx, ref, alt):
        """Extract AF and AD for Strelka from cyvcf2 variant."""
        try:
            ref_u = f"{ref}U"
            alt_u = f"{alt.split(',')[0]}U"  # Take first ALT allele

            ref_counts = (
                variant.format(ref_u)[sample_idx] if ref_u in variant.FORMAT else None
            )
            alt_counts = (
                variant.format(alt_u)[sample_idx] if alt_u in variant.FORMAT else None
            )

            if ref_counts is not None and alt_counts is not None:
                # Strelka uses tier1,tier2 format
                tier1_ref = int(ref_counts[0]) if len(ref_counts) > 0 else 0
                tier1_alt = int(alt_counts[1]) if len(alt_counts) > 1 else 0

                af = (
                    float(tier1_alt) / (tier1_alt + tier1_ref)
                    if (tier1_alt + tier1_ref) > 0
                    else 0.0
                )

                # Calculate AD as sum of all nucleotide counts (tier1 only)
                au = (
                    int(variant.format("AU")[sample_idx][0])
                    if "AU" in variant.FORMAT
                    else 0
                )
                cu = (
                    int(variant.format("CU")[sample_idx][0])
                    if "CU" in variant.FORMAT
                    else 0
                )
                gu = (
                    int(variant.format("GU")[sample_idx][0])
                    if "GU" in variant.FORMAT
                    else 0
                )
                tu = (
                    int(variant.format("TU")[sample_idx][0])
                    if "TU" in variant.FORMAT
                    else 0
                )
                ad = int(au + cu + gu + tu)

                return af, ad
        except (IndexError, TypeError):
            pass

        return None, None

    @staticmethod
    def _extract_somaticsniper_af_ad_cyvcf2(variant, sample_idx):
        """Extract AF and AD for SomaticSniper from cyvcf2 variant."""
        try:
            dp4 = variant.format("DP4")[sample_idx] if "DP4" in variant.FORMAT else None

            if dp4 is not None and len(dp4) == 4:
                ref_counts = int(dp4[0]) + int(dp4[1])
                alt_counts = int(dp4[2]) + int(dp4[3])

                total = ref_counts + alt_counts
                af = float(alt_counts) / total if total > 0 else 0.0

                return af, int(alt_counts)
        except (IndexError, TypeError):
            pass

        return None, None

    @staticmethod
    def create_vcf_index(vcf_path: str) -> None:
        """
        Create index for vcf file using tabix.
        """
        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"File not found: {vcf_path}")

        if not os.path.exists(f"{vcf_path}.tbi"):
            run(["tabix", "-p", "vcf", vcf_path], check=True)
        else:
            print("Index already exists.")

    @staticmethod
    def extract_sample_names_from_vcf(vcf_path: str) -> list:
        """
        Extracts sample names from vcf file using cyvcf2.
        """
        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"File not found: {vcf_path}")

        vcf = VCF(vcf_path)
        sample_names = vcf.samples
        vcf.close()

        return list(sample_names) if sample_names else []

    @staticmethod
    def extract_variantcaller_from_vcf(vcf_path: str) -> Optional[str]:
        """
        Extracts the variant caller from the VCF file using multiple detection methods
        with the pysam library for robust and efficient parsing.

        Args:
            vcf_path: The path to the VCF file.

        Returns:
            The name of the variant caller as a string (e.g., "mutect", "varscan")
            or None if no caller can be identified.
        """
        # List of known caller identifiers
        KNOWN_CALLERS = [
            "mutect",
            "varscan",
            "vardict",
            "haplotypecaller",
            "strelka",
            "varnet",
            "freebayes",
        ]

        try:
            with pysam.VariantFile(vcf_path) as vcf_file:
                # Method 1: Check '##source' header lines
                for record in vcf_file.header.records:
                    if record.key == "source":
                        source_value = str(record.value).lower()

                        for caller in KNOWN_CALLERS:
                            if caller in source_value:
                                return caller
                        break

                # Method 2: Check command line information in header
                for record in vcf_file.header.records:
                    if hasattr(record, "key") and record.key in [
                        "command",
                        "commandline",
                    ]:
                        command_value = str(record.value).lower()
                        for caller in KNOWN_CALLERS:
                            if caller in command_value:
                                return caller

                # Method 3: Check other header records for caller names
                header_text = str(vcf_file.header).lower()
                for caller in KNOWN_CALLERS:
                    if caller in header_text:
                        return caller

                # Method 4: Use signatures specifically for haplotypecaller vs mutect distinction
                info_fields = set(
                    field.lower() for field in vcf_file.header.info.keys()
                )
                format_fields = set(
                    field.lower() for field in vcf_file.header.formats.keys()
                )
                filter_fields = set(
                    field.lower() for field in vcf_file.header.filters.keys()
                )

                # Check for definitive signatures to distinguish haplotypecaller vs mutect
                if "tlod" in info_fields:  # Tumor LOD - unique to Mutect2
                    return "mutect"
                elif "germq" in info_fields:  # Germline quality - Mutect2
                    return "mutect"
                elif (
                    "strand_bias" in filter_fields and "weak_evidence" in filter_fields
                ):  # Mutect2 filters
                    return "mutect"
                elif (
                    "min_dp" in format_fields
                ):  # GVCF field - unique to HaplotypeCaller
                    return "haplotypecaller"
                elif (
                    "rgq" in format_fields
                ):  # Reference genotype quality - HaplotypeCaller
                    return "haplotypecaller"
                elif "failed_hc_filter" in filter_fields:  # Specific HC filter
                    return "haplotypecaller"

        except (IOError, ValueError, AttributeError) as e:
            logger.error(f"Error reading VCF file {vcf_path}: {e}")
            return None

        return None
