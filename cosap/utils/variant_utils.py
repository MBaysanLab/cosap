import tempfile
from pathlib import Path

import pandas as pd


class VariantUtils:
    @staticmethod
    def convert_list_to_ensembl_vep_input(variants: list, output_dir: str) -> str:
        """
        Converts list of variants to default vep input format and writes to a temporary file.
        The default vep input format is:
            1   881907    881906    -/C   +
            2   946507    946507    G/C   +
        """

        try:
            # Create directory using pathlib (more modern)
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            # Create temporary file
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".txt", dir=output_dir
            ) as f:
                for variant in variants:
                    end = variant["POS"] + len(variant["REF"]) - 1
                    f.write(
                        f"{variant['CHROM']}\t{variant['POS']}\t{end}\t{variant['REF']}/{variant['ALT']}\t+\n"
                    )
            return f.name
        except PermissionError:
            print(f"Permission denied to create directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory or file: {e}")

    @staticmethod
    def convert_list_to_vcf_input(variants: list, output_dir: str) -> str:
        """
        Converts list of variants to default vcf input format and writes to a temporary file.
        The default vcf input format is:
            ##fileformat=VCFv4.2
            #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
            1       881907  .       -       C       .       .       .
            2       946507  .       G       C       .       .       .

        Skips variants with <NON_REF> ALT alleles.
        Only includes variants from main GRCh38 contigs (chr1-22, X, Y, M).
        """
        
        # Define main GRCh38 contigs (chromosomes 1-22, X, Y, MT/M)
        MAIN_CONTIGS = {
            '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
            '21', '22', 'X', 'Y', 'M', 'MT',
            'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
            'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
            'chr21', 'chr22', 'chrX', 'chrY', 'chrM', 'chrMT'
        }

        try:
            # Create directory using pathlib (more modern)
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            # Create temporary file
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".vcf", dir=output_dir
            ) as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                
                skipped_non_ref = 0
                skipped_contigs = 0
                written_variants = 0
                
                for variant in variants:
                    # Skip variants with <NON_REF> ALT alleles
                    if variant['ALT'] == '<NON_REF>':
                        skipped_non_ref += 1
                        continue
                    
                    # Skip variants not in main contigs (decoy, alt, patches, etc.)
                    if variant['CHROM'] not in MAIN_CONTIGS:
                        skipped_contigs += 1
                        continue
                    
                    variant_id = f"{variant['CHROM']}_{variant['POS']}_{variant['REF']}_{variant['ALT']}"
                    f.write(
                        f"{variant['CHROM']}\t{variant['POS']}\t{variant_id}\t{variant['REF']}\t{variant['ALT']}\t.\t.\t.\n"
                    )
                    written_variants += 1
                
                # Log statistics
                if skipped_non_ref > 0:
                    print(f"Skipped {skipped_non_ref} variants with <NON_REF> ALT alleles")
                if skipped_contigs > 0:
                    print(f"Skipped {skipped_contigs} variants from non-main contigs (decoy/alt/patches)")
                print(f"Written {written_variants} variants to {f.name}")
                
            return f.name
        except PermissionError:
            print(f"Permission denied to create directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory or file: {e}")

    @staticmethod
    def convert_list_to_annovar_input(variants: list, output_dir: str) -> str:
        """
        Converts list of variants to default annovar input format and writes to a temporary file.
        The default annovar input format is:
            1 948921 948921 T C
            1 1404001 1404001 G T
    
        Skips variants with <NON_REF> ALT alleles.
        """
        try:
            # Create directory using pathlib (more modern)
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            # Create temporary file
            with tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".txt", dir=output_dir
            ) as f:
                for variant in variants:
                    # Skip variants with <NON_REF> ALT alleles
                    if variant['ALT'] == '<NON_REF>':
                        continue
                    
                    chrom = variant["CHROM"].replace("chr", "")
                    end = variant["POS"] + len(variant["REF"]) - 1
                    f.write(
                        f"{chrom}\t{variant['POS']}\t{end}\t{variant['REF']}\t{variant['ALT']}\n"
                    )
            return f.name
        except PermissionError:
            print(f"Permission denied to create directory: {output_dir}")
        except OSError as e:
            print(f"Error creating directory or file: {e}")

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
