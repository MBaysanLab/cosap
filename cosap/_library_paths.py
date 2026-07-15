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
    return os.path.join(
        AppConfig.LIBRARY_PATH,
        *filename,
    )


def _get_file_path_hg19(*filename):
    return os.path.join(
        AppConfig.LIBRARY_PATH,
        "ref_genome_indexes",
        "hg19_bundle",
        *filename,
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
    GENEFUSE_CANCER_GENE_LIST: str = _get_file_path_hg38("genefuse_cancer.hg38.csv")
    MSISENSOR_MICROSATELLITES: str = _get_file_path_hg38("msisensor_hg38.list")
    CNVKIT_ANNOTATION: str = _get_file_path_hg38("cnvkit_hg38_refFlat.txt")
    CNVKIT_ACCESS: str = _get_file_path_hg38("cnvkit_hg38_access-5kb.bed")
    ANNOTSV: str = _get_file_path_hg38("annotsv", "bin", "AnnotSV")
    CLASSIFYCNV: str = _get_file_path_hg38("classifycnv")
    EXOMISER_JAR: str = _get_file_path_hg38(
        "exomiser-cli-14.1.0", "exomiser-cli-14.1.0.jar"
    )
    METARNN = _get_file_path_hg38("MetaRNN")


@dataclass
class _LibraryPaths38(_LibraryPaths):
    REF_DIR: str = AppConfig.LIBRARY_PATH
    REF_FASTA: str = _get_file_path_hg38("Homo_sapiens_assembly38.fasta")
    REF_DICT: str = _get_file_path_hg38("Homo_sapiens_assembly38.dict")
    REF_ELFASTA: str = _get_file_path_hg38("Homo_sapiens_assembly38.elfasta")
    DBSNP: str = _get_file_path_hg38("Homo_sapiens_assembly38.dbsnp138.vcf")
    DBSNP_ELSITES: str = _get_file_path_hg38("Homo_sapiens_assembly38.dbsnp138.elsites")
    MILLS_INDEL: str = _get_file_path_hg38(
        "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    )
    MILLS_INDEL_ELSITES: str = _get_file_path_hg38(
        "Mills_and_1000G_gold_standard.indels.hg38.elsites"
    )
    COSMIC: str = _get_file_path_hg38("cosmic_hg19_lifted_over.vcf")
    ONE_THOUSAND_G: str = _get_file_path_hg38(
        "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    )
    ONE_THOUSAND_G_ELSITES: str = _get_file_path_hg38(
        "1000G_phase1.snps.high_confidence.hg38.elsites"
    )
    BWA_ASSEMBLY: str = _get_file_path_hg38("Homo_sapiens_assembly38.fasta")
    BOWTIE2_ASSEMBLY: str = _get_file_path_hg38("Homo_sapiens_assembly38")
    CHAIN_FILE: str = _get_file_path_hg38("hg19ToHg38.over.chain.gz")
    # GermlineCNVCaller reference files
    CONTIG_PLOIDY_PRIORS: str = _get_file_path_hg38("contig_ploidy_priors_hg38.tsv")
    # GATK mitochondrial Mutect2 resource bundle.
    MT_DIR: str = _get_file_path_hg38("chrM")
    MT_FASTA: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta")
    MT_DICT: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.dict")
    MT_FASTA_INDEX: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta.fai")
    MT_FASTA_AMB: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta.amb")
    MT_FASTA_ANN: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta.ann")
    MT_FASTA_BWT: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta.bwt")
    MT_FASTA_PAC: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta.pac")
    MT_FASTA_SA: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.fasta.sa")
    MT_SHIFTED_FASTA: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta")
    MT_SHIFTED_DICT: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict")
    MT_SHIFTED_FASTA_INDEX: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai")
    MT_SHIFTED_FASTA_AMB: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb")
    MT_SHIFTED_FASTA_ANN: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann")
    MT_SHIFTED_FASTA_BWT: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt")
    MT_SHIFTED_FASTA_PAC: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac")
    MT_SHIFTED_FASTA_SA: str = _get_file_path_hg38("chrM", "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa")
    MT_SHIFT_BACK_CHAIN: str = _get_file_path_hg38("chrM", "ShiftBack.chain")
    MT_BLACKLIST: str = _get_file_path_hg38("chrM", "blacklist_sites.hg38.chrM.bed")
    MT_BLACKLIST_INDEX: str = _get_file_path_hg38("chrM", "blacklist_sites.hg38.chrM.bed.idx")


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
