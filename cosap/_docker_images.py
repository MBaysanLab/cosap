from dataclasses import dataclass


@dataclass
class DockerImages:
    DEEPVARIANT: str = "google/deepvariant"
    ENSEMBL_VEP: str = "ensemblorg/ensembl-vep"
    PARABRICKS: str = "nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1"
    COSAP: str = "itubioinformatics/cosap"
    PHARMCAT: str = "pgkb/pharmcat"
    STRELKA2: str = "iarcbioinfo/strelka2-nf"
    MANTA: str = "szarate/manta:v1.6.0"
    VARNET: str = "kiranchari/varnet"
