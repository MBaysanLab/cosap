from dataclasses import dataclass

@dataclass
class DockerContainers:
    DEEPVARIANT: str = "google/deepvariant"
    COSAP: str = "itubioinformatics/cosap"