from cosap.utils.path_utils import PathUtils
from cosap.utils.vcf_utils import VCFUtils
from cosap.utils.variant_utils import VariantUtils

import gzip
import os
from subprocess import run

import pandas as pd
import pyranges as pr
from pandas.api.types import is_numeric_dtype


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


def convert_vcf_to_json(path: str, caller_type: str = "mutect", sample_name: str = "TUMOR") -> list:
    return VCFUtils.convert_vcf_to_json(path, caller_type, sample_name)


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
