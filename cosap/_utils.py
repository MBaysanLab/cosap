import gzip
import os
import psutil

import pandas as pd
from pandas.api.types import is_numeric_dtype
from subprocess import run
import gzip

from cosap.utils.path_utils import PathUtils
from cosap.utils.variant_utils import VariantUtils
from cosap.utils.vcf_utils import VCFUtils


def join_paths(path: str, *paths) -> str:
    """Joins and normalizes paths according to os standards"""
    return PathUtils.join_paths(path, *paths)


def convert_to_absolute_path(path: str) -> str:
    """Converts path to absolute path"""
    return PathUtils.convert_to_absolute_path(path)


def read_vcf_into_df(path: str) -> pd.DataFrame:
    return VCFUtils.read_vcf_into_df(path)


def convert_vcf_to_tsv(path: str, caller_type: str = "mutect") -> str:
    return VCFUtils.convert_vcf_to_tsv(path, caller_type)


def calculate_strelka_af(row: pd.Series) -> float:
    return VCFUtils.calculate_strelka_af(row)


def calculate_strelka_ad(row: pd.Series) -> float:
    return VCFUtils.calculate_strelka_ad(row)


def rename_columns_strelka(vcf_df, sample_name):
    return VCFUtils.get_column_mappings_strelka(vcf_df, sample_name)


def rename_columns_haplotypecaller(vcf_df, sample_name):
    return VCFUtils.get_column_mappings_haplotypecaller(vcf_df, sample_name)


def rename_columns_varscan(vcf_df, sample_name):
    return VCFUtils.get_column_mappings_varscan(vcf_df, sample_name)


def rename_columns_generic(vcf_df, sample_name):
    return VCFUtils.get_column_mappings_generic(vcf_df, sample_name)


def rename_columns_varnet(vcf_df, sample_name):
    return VCFUtils.get_column_mappings_varnet(vcf_df, sample_name)



def convert_vcf_to_json(
    path: str, caller_type: str = "mutect", sample_name: str = "TUMOR", filter: str = "PASS"
) -> list:
    return VCFUtils.convert_vcf_to_json(path, caller_type, sample_name, filter)


def is_valid_path(path: str) -> bool:
    """
    Returns True if path exists and is not empty.
    """
    return PathUtils.is_valid_path(path)


def get_commonpath_from_config(config: dict) -> str:
    """
    Returns the commonpath of paths that are in the config.
    """
    return PathUtils.get_commonpath_from_config(config)


def convert_list_to_ensembl_vep_input(variants: list, workdir: str) -> str:
    return VariantUtils.convert_list_to_ensembl_vep_input(variants, workdir)


def convert_list_to_annovar_input(variants: list, workdir: str) -> str:
    return VariantUtils.convert_list_to_annovar_input(variants, workdir)


def get_variants_within_bed_regions(variants_df, bed_df):
    return VariantUtils.get_variants_within_bed_regions(variants_df, bed_df)

class MultipleFilesFoundError(Exception):
    """Exception raised when multiple files are found when only one is expected."""
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


def get_bam_index_path(bam_path: str) -> str:
    root_path, extension = os.path.splitext(bam_path)

    if extension != ".bam":
        raise ValueError(f"Extension of BAM file is '{extension}'.")

    if not os.path.isfile(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    possible_index_files = [
        root_path + "." + "bam" + "." + "bai",
        root_path + "." + "bai",
    ]

    present_index_files = [
        file for file in possible_index_files
        if os.path.isfile(file)
    ]

    if len(present_index_files) == 0:
        raise FileNotFoundError(f"No BAM index file found for: {bam_path}")
    
    if len(present_index_files) > 1:
        raise MultipleFilesFoundError(f"Multiple possible BAM index files found for: {bam_path}")
    
    bai_path = present_index_files[0]

    return bai_path


def swap_dict_keys_and_values(a_dict: dict) -> dict:
    return {value: key for key, value in a_dict.items()}


def prompt_continue(message: str = "Continue?"):
    while True:
        user_input = input(f"{message} [Y/n] ").strip().lower()
        if user_input in {"y", "yes", ""}:
            return True
        if user_input in {"n", "no"}:
            return False
        print("Invalid input. Please enter 'yes' or 'no'.")


def current_available_system_memory():
    current_system_memory_info = psutil.virtual_memory()
    current_available_system_memory = current_system_memory_info.available
    return current_available_system_memory