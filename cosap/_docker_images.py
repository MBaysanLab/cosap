from dataclasses import dataclass


@dataclass
class DockerImages:
    DEEPVARIANT: str = "google/deepvariant:1.10.0-beta"
    ENSEMBL_VEP: str = "ensemblorg/ensembl-vep:release_115.2"
    PARABRICKS: str = "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"
    COSAP: str = "itubioinformatics/cosap"
    PHARMCAT: str = "pgkb/pharmcat"
    STRELKA2: str = "iarcbioinfo/strelka2-nf"
    MANTA: str = "szarate/manta:v1.6.0"
    GATK: str = "broadinstitute/gatk"
    HAPLOCHECKER: str = "us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124"
    VARNET: str = "kiranchari/varnet"
    METARNN: str = "itubioinformatics/metarnn"
