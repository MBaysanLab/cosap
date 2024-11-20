import gzip
import os
from subprocess import run

import pandas as pd
from pandas.api.types import is_numeric_dtype


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
            awk_command = rf"zcat {path} | awk -F '\t' '/^#/ || $4 != $5' > {new_path}"
        else:
            awk_command = rf"awk -F '\t' '/^#/ || $4 != $5' {path} > "

        run(["bash", "-c", awk_command], check=True)
        return new_path

    @staticmethod
    def convert_vcf_to_tsv(path: str, caller_type: str = "mutect") -> str:
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
        elif caller_type == "haplotypecaller":
            return {
                "AD": f"{sample_name}.AD",
                "DP": f"{sample_name}.DP",
            }
        elif caller_type == "varscan":
            vcf_df["varscan_AF"] = (
                vcf_df[f"{sample_name}.FREQ"].str.replace("%", "").astype(float) / 100
            )

            return {
                "AF": "varscan_AF",
            }
        elif caller_type == "mutect2":
            return {
                "AF": f"{sample_name}.AF",
                "AD": f"{sample_name}.AD",
                "DP": f"{sample_name}.DP",
            }
        elif caller_type == "varnet":
            return {
                "AF": "SAMPLE.AF",
                "AD": "SAMPLE.AO",
                "DP": "SAMPLE.DP",
            }
        else:
            raise ValueError(f"Unknown tool type: {caller_type}")

    @staticmethod
    def convert_vcf_to_json(
        path: str, caller_type: str = "mutect2", sample_name: str = "TUMOR"
    ) -> list:
        """
        Returns list of variants as json objects.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

        variants_table = VCFUtils.convert_vcf_to_tsv(path, caller_type)
        vcf_df = pd.read_csv(variants_table, sep="\t")

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

        if not is_numeric_dtype(final_df["AD"]):
            final_df["AD"] = final_df["AD"].apply(
                lambda x: x.split(",")[-1] if "," in x else x
            )

        return final_df.to_dict("records")
