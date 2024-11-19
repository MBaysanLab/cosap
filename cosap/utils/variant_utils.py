
import tempfile
import pandas as pd

class VariantUtils:
    @staticmethod
    def convert_list_to_ensembl_vep_input(variants: list, workdir: str) -> str:
        """
        Converts list of variants to default vep input format and writes to a temporary file.
        The default vep input format is:
            1   881907    881906    -/C   +
            2   946507    946507    G/C   +
        """
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".txt", dir=workdir
        ) as f:
            for variant in variants:

                end = variant["POS"] + len(variant["REF"]) - 1

                f.write(
                    f"{variant['CHROM']}\t{variant['POS']}\t{end}\t{variant['REF']}/{variant['ALT']}\t+\n"
                )
            return f.name

    @staticmethod
    def convert_list_to_annovar_input(variants: list, workdir: str) -> str:
        """
        Converts list of variants to default annovar input format and writes to a temporary file.
        The default annovar input format is:
            1 948921 948921 T C
            1 1404001 1404001 G T
        """
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".avinput", dir=workdir
        ) as f:
            for variant in variants:
                chrom = variant["CHROM"].replace("chr", "")
                end = variant["POS"] + len(variant["REF"]) - 1
                f.write(
                    f"{chrom}\t{variant['POS']}\t{end}\t{variant['REF']}\t{variant['ALT']}\n"
                )
            return f.name

    @staticmethod
    def get_variants_within_bed_regions(variants_df, bed_df):
        """
        Returns variants that are within the bed regions.
        """
        variants_df["POS"] = variants_df["POS"].astype(int)
        bed_df["POS"] = bed_df["POS"].astype(int)
        bed_df["END"] = bed_df["END"].astype(int)
        variants_df = pd.merge(
            variants_df,
            bed_df,
            on=["CHROM"],
            how="inner",
            suffixes=("_variants", "_bed"),
        )
        variants_df = variants_df[
            (variants_df["POS_variants"] >= variants_df["POS_bed"])
            & (variants_df["POS_variants"] <= variants_df["END"])
        ]
        return variants_df