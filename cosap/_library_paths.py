from __future__ import annotations

import os
from dataclasses import dataclass
from threading import Lock

from ._config import AppConfig


class _LibraryMeta(type):
    _instances = {}
    _lock: Lock = Lock()

    def __call__(cls, grch: str = "hg38") -> _LibraryPaths:
        with cls._lock:
            if cls not in cls._instances:
                if grch == "hg38":
                    instance = _LibraryPaths38
                else:
                    instance = _LibraryPaths19
                cls._instances[cls] = instance.__call__()
        return cls._instances[cls]


@dataclass
class _LibraryPaths:
    # PICARD: str = os.path.join(AppConfig.LIBRARY_PATH, "picard.jar")
    # GATK: str = os.path.join(AppConfig.LIBRARY_PATH, "GenomeAnalysisTK.jar")
    # GATK4: str = os.path.join(AppConfig.LIBRARY_PATH, "gatk-4.1.0.0", "gatk")
    # VARSCAN: str = os.path.join(AppConfig.LIBRARY_PATH, "VarScan.v2.3.9.jar")
    # FASTQC: str = os.path.join(AppConfig.LIBRARY_PATH, "FastQC", "fastqc")
    # FASTP: str = os.path.join(AppConfig.LIBRARY_PATH, "fastp", "fastp")
    # STRELKA: str = os.path.join(
    #     AppConfig.LIBRARY_PATH,
    #     "strelka-2.9.10.centos6_x86_64",
    #     "bin",
    #     "configureStrelkaSomaticWorkflow.py",
    # )
    # SOMATICSNIPER: str = os.path.join(
    #     AppConfig.LIBRARY_PATH, "somatic-sniper", "build", "bin", "bam-somaticsniper"
    # )


@dataclass
class _LibraryPaths38(_LibraryPaths):
    REF_DIR: str = os.path.join(
        AppConfig.LIBRARY_PATH, "ref_genome_indexes", "hg38_bundle"
    )
    DBSNP: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg38_bundle",
        "dbsnp_146.hg38.vcf.gz",
    )
    MILLS_INDEL: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg38_bundle",
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    )
    COSMIC: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg19_bundle",
        "cosmic_hg19_lifted_over.vcf",
    )
    ANNOVAR: str = os.path.join(AppConfig.LIBRARY_PATH, "annovar")
    ANNOVAR_DB: str = os.path.join(AppConfig.LIBRARY_PATH, "annovar", "humandb_38")
    ONE_THOUSAND_G: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg38_bundle",
        "1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    )


@dataclass
class _LibraryPaths19(_LibraryPaths):
    REF_DIR: str = os.path.join(
        AppConfig.LIBRARY_PATH, "ref_genome_indexes", "hg19_bundle"
    )
    DBSNP: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg19_bundle",
        "dbsnp_138.hg19.vcf.gz",
    )
    MILLS_INDEL: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg19_bundle",
        "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
    )
    COSMIC: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg19_bundle",
        "cosmic_hg19_lifted_over.vcf",
    )
    ANNOVAR: str = os.path.join(AppConfig.LIBRARY_PATH, "annovar")
    ANNOVAR_DB: str = os.path.join(AppConfig.LIBRARY_PATH, "annovar", "humandb")
    ONE_THOUSAND_G: str = os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg19_bundle",
        "1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz",
    )


@dataclass
class LibraryPaths(metaclass=_LibraryMeta):
    """
    Only this class is imported from outside
    Metaclass will automatically create the correct instance
    """