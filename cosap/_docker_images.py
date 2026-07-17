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
    FREEC: str = (
        "quay.io/biocontainers/control-freec:11.6b--hde5307d_3@sha256:488fa257f3aa201bba1d8e34b88ea499dde33ac2143df023d7046ea768fbf126"
    )
    CLINCNV: str = "itubioinformatics/clincnv:latest"
    HAPLOCHECKER: str = "us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124"
    VARNET: str = "kiranchari/varnet"
    METARNN: str = "itubioinformatics/metarnn"
