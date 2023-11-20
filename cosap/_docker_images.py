from dataclasses import dataclass


@dataclass
class DockerImages:
    DEEPVARIANT: str = "google/deepvariant"
    ENSEMBL_VEP: str = "ensemblorg/ensembl-vep"
    PARABRICKS: str = "nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1"
    COSAP: str = "itubioinformatics/cosap"
    PHARMCAT: str = "pgkb/pharmcat"
