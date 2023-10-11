import os

import pandas as pd


def join_paths(path: str, *paths) -> str:
    """Joins and normalizes paths according to os standards"""
    return os.path.normpath(os.path.join(path, *paths))


def read_vcf_into_df(path: str) -> pd.DataFrame:
    import io

    with open(path, "r") as f:
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
        ).rename(
            columns={"#Chr": "Chr", "Ref.Gene": "Gene", "Func.refGene": "Function"}
        )

    df.reset_index(inplace=True)
    df.rename(columns={"index": "id"}, inplace=True)
    return df


def convert_vcf_to_json(path: str) -> list:
    """
    Returns list of variants as json objects.
    """
    vcf_df = read_vcf_into_df(path)
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
