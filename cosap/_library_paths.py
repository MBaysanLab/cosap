from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from threading import Lock

from ._config import AppConfig


class _LibraryMeta(type):
    _instances = {}
    _lock: Lock = Lock()

    def __call__(cls) -> _LibraryPaths:

        version = AppConfig.REF_VERSION
        with cls._lock:
            if cls not in cls._instances:
                if version.lower() == "hg38":
                    instance = _LibraryPaths38
                elif version.lower() == "hg19":
                    instance = _LibraryPaths19
                cls._instances[cls] = instance.__call__()
        return cls._instances[cls]


def _get_file_path_hg38(*filename):
    return (
        os.path.join(
            AppConfig.LIBRARY_PATH,
            *filename,
        )
    )


def _get_file_path_hg19(*filename):
    return (
        os.path.join(
            AppConfig.LIBRARY_PATH,
            "ref_genome_indexes",
            "hg19_bundle",
            *filename,
        )
    )


@dataclass
class _LibraryPaths:
    ENSEMBL_VEP: str = _get_file_path_hg38("ensembl-vep")
    ANNOVAR: str = _get_file_path_hg38("annovar")
    INTERVAR: str = _get_file_path_hg38("intervar")
    CANCERVAR: str = _get_file_path_hg38("cancervar")
    PHARMCAT_DIR: str = _get_file_path_hg38("pharmcat")
    PHARMCAT_PREPROCESSOR: str = os.path.join(
        PHARMCAT_DIR, "preprocessor", "pharmcat_vcf_preprocessor.py"
    )
    PHARMCAT_JAR: str = os.path.join(PHARMCAT_DIR, "pharmcat-2.8.3-all.jar")
    VARNET: str = _get_file_path_hg38("varnet")
    GENEFUSE_CANCER_GENE_LIST: str = _get_file_path_hg38("genefuse_cancer.hg38.csv")
    MSISENSOR_MICROSATELLITES: str = _get_file_path_hg38("msisensor_hg38.list")
    CNVKIT_ANNOTATION: str = _get_file_path_hg38("cnvkit_hg38_refFlat.txt")
    CNVKIT_ACCESS: str = _get_file_path_hg38("cnvkit_hg38_access-5kb.bed")
    ANNOTSV: str = _get_file_path_hg38("annotsv", "bin", "AnnotSV")
    CLASSIFYCNV: str = _get_file_path_hg38("classifycnv")


@dataclass
class _LibraryPaths38(_LibraryPaths):
    REF_DIR: str = AppConfig.LIBRARY_PATH
    REF_FASTA: str = _get_file_path_hg38(
        "Homo_sapiens_assembly38.fasta"
    )
    REF_ELFASTA: str = _get_file_path_hg38(
        "Homo_sapiens_assembly38.elfasta"
    )
    DBSNP: str = _get_file_path_hg38(
        "Homo_sapiens_assembly38.dbsnp138.vcf"
    )
    DBSNP_ELSITES: str = _get_file_path_hg38(
        "Homo_sapiens_assembly38.dbsnp138.elsites"
    )
    MILLS_INDEL: str = _get_file_path_hg38(
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    )
    MILLS_INDEL_ELSITES: str = _get_file_path_hg38(
        "Mills_and_1000G_gold_standard.indels.hg38.elsites"
    )
    COSMIC: str = _get_file_path_hg38(
        "cosmic_hg19_lifted_over.vcf"
    )
    ONE_THOUSAND_G: str = _get_file_path_hg38(
        "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    )
    ONE_THOUSAND_G_ELSITES: str = _get_file_path_hg38(
        "1000G_phase1.snps.high_confidence.hg38.elsites"
    )
    BWA_ASSEMBLY: str = _get_file_path_hg38("Homo_sapiens_assembly38.fasta")
    BOWTIE2_ASSEMBLY: str = _get_file_path_hg38("Homo_sapiens_assembly38")


@dataclass
class _LibraryPaths19(_LibraryPaths):
    REF_DIR: str = _get_file_path_hg19()
    DBSNP: str = _get_file_path_hg19(
        "dbsnp_138.hg19.vcf.gz",
    )
    MILLS_INDEL: str = _get_file_path_hg19(
        "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
    )
    COSMIC: str = _get_file_path_hg19(
        "cosmic_hg19_lifted_over.vcf",
    )
    ONE_THOUSAND_G: str = _get_file_path_hg19(
        "1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz",
    )


@dataclass
class LibraryPaths(metaclass=_LibraryMeta):
    """
    Only this class is imported from outside
    Metaclass will automatically create the correct instance
    """

    # Check if all the paths are correct
    def __post_init__(self):
        for key, value in self.__dict__.items():
            if not os.path.exists(value):
                warnings.warn(f"Path {value} does not exist")
