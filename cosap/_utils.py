import os

import pandas as pd
from pandas.api.types import is_numeric_dtype
from subprocess import run
import pyranges as pr
import gzip


def join_paths(path: str, *paths) -> str:
    """Joins and normalizes paths according to os standards"""
    return os.path.normpath(os.path.join(path, *paths))


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


def convert_vcf_to_tsv(path: str, caller_type: str = "mutect") -> str:

    # Gatk VariantsToTable gives error when ref and alt alleles are same, so remove those lines
    awk_command = rf"awk -F '\t' '/^#/ || $4 != $5' {path} > {path}.tmp && mv {path}.tmp {path}"
    run(["bash", "-c", awk_command], check=True)

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
        "INFO",
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

def calculate_strelka_ad(row: pd.Series) -> float:
    """
    Calculates AD as sum of AU, CU, GU, TU columns.
    """
    return sum([int(i) for i in row[["TUMOR.AU", "TUMOR.CU", "TUMOR.GU", "TUMOR.TU"]].values[0].split(",")])

def convert_vcf_to_json(
    path: str, caller_type: str = "mutect", sample_name: str = "TUMOR"
) -> list:
    """
    Returns list of variants as json objects.
    """

    sample_name = sample_name.upper()

    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    variants_table = convert_vcf_to_tsv(path, caller_type)

    vcf_df = pd.read_csv(variants_table, sep="\t")

    if caller_type == "strelka":
        for i, row in vcf_df.iterrows():
            af = calculate_strelka_af(row)
            vcf_df.at[i, "AF"] = af
            ad = calculate_strelka_ad(row)
            vcf_df.at[i, "AD"] = ad

    # Convert all columns to UPPERCASE
    vcf_df.columns = vcf_df.columns.str.upper()

    # Rename columns
    if f"{sample_name}.AF" in vcf_df.columns:
        vcf_df.rename(
            columns={
                f"{sample_name}.AF": "AF",
                f"{sample_name}.AD": "AD",
                f"{sample_name}.DP": "DP",
            },
            inplace=True,
        )
    if "SAMPLE.AF" in vcf_df.columns:
        vcf_df.rename(
            columns={
                "SAMPLE.AF": "AF",
                "SAMPLE.AD": "AD",
                "SAMPLE.DP": "DP",
            },
            inplace=True,
        )
    # VarScan uses FREQ instead of AF
    if f"{sample_name}.FREQ" in vcf_df.columns:
        vcf_df.rename(
            columns={
                f"{sample_name}.FREQ": "AF",
            },
            inplace=True,
        )
        vcf_df["AF"] = vcf_df["AF"].replace("%", "", regex=True).astype(float) / 100

    if caller_type.lower() == "varnet":
        vcf_df.rename(
            columns={
                f"SAMPLE.AO": "AD",
            },
            inplace=True,
        )

    columns_to_keep = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "AF",
        "AD",
        "DP",
        "STATUS",
    ]

    vcf_df = vcf_df[columns_to_keep]

    # AD column contains comma separated values of ref and alt allele counts, get only alt allele counts
    # IF AD is str type and contains comma separated values, get the second value
    if not is_numeric_dtype(vcf_df["AD"]):
        vcf_df["AD"] = vcf_df["AD"].apply(lambda x: x.split(",")[-1] if "," in x else x)

    

    return vcf_df.to_dict("records")

def is_valid_path(path: str) -> bool:
    """
    Returns True if path exists and is not empty.
    """
    return os.path.exists(path) or os.path.isabs(path)


def get_commonpath_from_config(config: dict) -> str:
    """
    Returns the commonpath of paths that are in the config.
    """
    paths = []
    for key, value in config.items():
        if isinstance(value, str):
            if is_valid_path(value):
                paths.append(value)
        elif isinstance(value, list):
            for i in value:
                if is_valid_path(i):
                    paths.append(i)
        elif isinstance(value, dict):
            for v in value.values():
                if is_valid_path(v):
                    paths.append(v)
        else:
            raise ValueError("Config value is not a string or list of strings.")

    return os.path.commonpath(paths)


def convert_list_to_ensembl_vep_input(variants: list, workdir: str) -> str:
    """
    Converts list of variants to default vep input format and writes to a temporary file.
    The default vep input format is:
        1   881907    881906    -/C   +
        2   946507    946507    G/C   +
    """
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt", dir=workdir) as f:
        for variant in variants:

            end = variant["POS"] + len(variant["REF"]) - 1

            f.write(
                f"{variant['CHROM']}\t{variant['POS']}\t{end}\t{variant['REF']}/{variant['ALT']}\t+\n"
            )
        return f.name


def convert_list_to_annovar_input(variants: list, workdir:str) -> str:
    """
    Converts list of variants to default annovar input format and writes to a temporary file.
    The default annovar input format is:
        1 948921 948921 T C
        1 1404001 1404001 G T
    """
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".avinput", dir=workdir) as f:
        for variant in variants:
            chrom = variant["CHROM"].replace("chr", "")
            end = variant["POS"] + len(variant["REF"]) - 1
            f.write(
                f"{chrom}\t{variant['POS']}\t{end}\t{variant['REF']}\t{variant['ALT']}\n"
            )
        return f.name

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
