from dataclasses import dataclass


@dataclass
class PipelineConfig:
    MAPPER_TYPE: str
    MAPPER_THREADS: str
    SAMPLE_TYPE: str
    FASTQ_TRIMMED: bool
    FASTQ_DIR: str
    PATIENT_ID: str

    VARIANT_CALLER_TYPE: str

    BAM_DIR: str
    VCF_OUTPUT_DIR: str

    USER_CONFIG: str

