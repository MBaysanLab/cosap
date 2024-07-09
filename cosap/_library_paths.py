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


@dataclass
class _LibraryPaths:
    ENSEMBL_VEP: str = os.path.join(AppConfig.LIBRARY_PATH, "ensembl-vep")
    ANNOVAR: str = os.path.join(AppConfig.LIBRARY_PATH, "annovar")
    INTERVAR: str = os.path.join(AppConfig.LIBRARY_PATH, "intervar")
    CANCERVAR: str = os.path.join(AppConfig.LIBRARY_PATH, "cancervar")
    PHARMCAT_DIR: str = os.path.join(AppConfig.LIBRARY_PATH, "pharmcat")
    PHARMCAT_PREPROCESSOR: str = os.path.join(
        PHARMCAT_DIR, "preprocessor", "pharmcat_vcf_preprocessor.py"
    )
    PHARMCAT_JAR: str = os.path.join(PHARMCAT_DIR, "pharmcat-2.8.3-all.jar")
    VARNET: str = os.path.join(AppConfig.LIBRARY_PATH, "varnet")
    ANNOTSV: str = os.path.join(AppConfig.LIBRARY_PATH, "annotsv", "bin", "AnnotSV")
    CLASSIFYCNV: str = os.path.join(AppConfig.LIBRARY_PATH, "classifycnv")


@dataclass
class _LibraryPaths38(_LibraryPaths):
    REF_DIR: str = AppConfig.LIBRARY_PATH
    REF_FASTA: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly38.fasta"
    )
    REF_ELFASTA: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly38.elfasta"
    )
    REF_BED: str = os.path.join(AppConfig.LIBRARY_PATH, "hg38_intervals.bed")
    REF_GFF3: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens.GRCh38.111.gff3.gz"
    )
    DBSNP: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "Homo_sapiens_assembly38.dbsnp138.vcf",
    )
    DBSNP_ELSITES = os.path.join(
        AppConfig.LIBRARY_PATH,
        "Homo_sapiens_assembly38.dbsnp138.elsites",
    )
    MILLS_INDEL: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    )
    MILLS_INDEL_ELSITES: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "Mills_and_1000G_gold_standard.indels.hg38.elsites",
    )
    COSMIC: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "cosmic_hg19_lifted_over.vcf",
    )
    ONE_THOUSAND_G: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    )
    ONE_THOUSAND_G_ELSITES: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "1000G_phase1.snps.high_confidence.hg38.elsites",
    )
    BWA_ASSEMBLY: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly38.fasta"
    )
    BOWTIE2_ASSEMBLY: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly38"
    )
    GENEFUSE_CANCER_GENE_LIST: str = os.path.join(
        AppConfig.LIBRARY_PATH, "genefuse_cancer.hg38.csv"
    )
    MSISENSOR_MICROSATELLITES: str = os.path.join(
        AppConfig.LIBRARY_PATH, "msisensor_hg38.list"
    )
    CNVKIT_ANNOTATION: str = os.path.join(
        AppConfig.LIBRARY_PATH, "cnvkit_hg38_refFlat.txt"
    )
    CNVKIT_ACCESS: str = os.path.join(
        AppConfig.LIBRARY_PATH, "cnvkit_hg38_access-5kb.bed"
    )
    INTERVALS: str = os.path.join(AppConfig.LIBRARY_PATH, "intervals")
    STRINGENCIES: str = os.path.join(AppConfig.LIBRARY_PATH, "stringencies")


@dataclass
class _LibraryPaths19(_LibraryPaths):
    REF_DIR: str = AppConfig.LIBRARY_PATH
    REF_FASTA: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly19.fasta"
    )
    REF_ELFASTA: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly19.elfasta"
    )
    DBSNP: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "Homo_sapiens_assembly19.dbsnp138.vcf",
    )
    MILLS_INDEL: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "Mills_and_1000G_gold_standard.indels.hg19.vcf.gz",
    )
    ONE_THOUSAND_G: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "1000G_phase1.snps.high_confidence.hg19.vcf.gz",
    )
    BWA_ASSEMBLY: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly19.fasta"
    )
    BOWTIE2_ASSEMBLY: str = os.path.join(
        AppConfig.LIBRARY_PATH, "Homo_sapiens_assembly19"
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
