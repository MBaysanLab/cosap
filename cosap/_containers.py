from dataclasses import dataclass


@dataclass
class DockerContainers:
    DEEPVARIANT: str = "google/deepvariant"
    ENSEMBL_VEP: str = "ensemblorg/ensembl-vep"
    COSAP: str = "itubioinformatics/cosap"
