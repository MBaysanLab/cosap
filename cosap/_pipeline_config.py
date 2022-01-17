from dataclasses import dataclass
from typing import Pattern


@dataclass
class PipelineKeys:
    WORKDIR: str = "workdir"
    LIBRARY_PATH: str = "library_path"
    CREATION_DATE: str = "creation-date"
    VERSION: str = "version"
    MAPPING: str = "mapping"
    SORTING: str = "sorting"
    INDEX: str = "index"
    MERGE: str = "merge"
    MDUP: str = "mdup"
    TRIM: str = "trim"
    CALIBRATE: str = "calibrate"
    ELPREP_PROCESS: str = "elprep"
    VARIANT_CALLING: str = "variant-calling"
    ANNOTATION: str = "annotation"
    FINAL_OUTPUT: str = "final-output"


@dataclass
class PipelineBaseKeys:
    LIBRARY: str = "library"
    PARAMS: str = "params"
    INPUT: str = "input"
    OUTPUT: str = "output"
    SNAKEMAKE_OUTPUT: str = "snakemake_output"


@dataclass
class MappingKeys(PipelineBaseKeys):
    SAMPLE_TYPE: str = "sample-type"
    READ_GROUP: str = "read_groups"
    RG_ID: str = "ID"
    RG_SM: str = "SM"
    RG_LB: str = "LB"
    RG_PL: str = "PL"
    RG_PU: str = "PU"


@dataclass
class TrimmingKeys(PipelineBaseKeys):
    pass


@dataclass
class SortingKeys(PipelineBaseKeys):
    INPUT: str = "input"
    OUTPUT: str = "output"
    BAM_DIR: str = "bam-dir"


@dataclass
class IndexingKeys(PipelineBaseKeys):
    INPUT: str = "input"
    OUTPUT: str = "output"
    BAM_DIR: str = "bam-dir"


@dataclass
class SplitKeys(PipelineBaseKeys):
    INPUT: str = "inputs"
    OUTPUT: str = "outputs"
    OUTPUT_DIR: str = "output-dir"


@dataclass
class MergingKeys(PipelineBaseKeys):
    INPUTS: str = "inputs"
    OUTPUT: str = "output"


@dataclass
class MDUPKeys(PipelineBaseKeys):
    INPUT: str = "input"
    OUTPUT: str = "output"
    METRICS: str = "metrics"
    OUTPUT_DIR: str = "output-dir"


@dataclass
class BaseRecalibratorKeys(PipelineBaseKeys):
    INPUT: str = "input"
    OUTPUT: str = "output"
    TABLE: str = "table"


@dataclass
class VariantCallingKeys(PipelineBaseKeys):
    GERMLINE_INPUT: str = "normal_input"
    TUMOR_INPUT: str = "tumor_input"
    GERMLINE_SAMPLE_NAME: str = "germline_sample_name"
    TUMOR_SAMPLE_NAME: str = "tumor_sample_name"
    UNFILTERED_VARIANTS_OUTPUT: str = "unfiltered_variants"
    SNP_OUTPUT: str = "snp_output"
    INDEL_OUTPUT: str = "indel_output"
    OTHER_VARIANTS_OUTPUT: str = "other_variants_output"
    PILEUPS: str = "pileups"
    CONDA_ENV: str = "conda_env"


@dataclass
class AnnotatorKeys(PipelineBaseKeys):
    pass

@dataclass
class ElprepKeys(MDUPKeys, BaseRecalibratorKeys):
    pass


@dataclass
class DefaultValues:
    DEFAULT_ENV: str = "default_environment"
    DEFAULT_ENV_PY2: str = "py2_environment"
    DEFAULT_GERMLINE_SAMPLE_NAME: str = "normal_sample"
    DEFAULT_TUMOR_SAMPLE_NAME: str = "tumor_sample"

@dataclass
class Py2ModuleConstraints:
    VARIANT_CALLERS: list[Pattern] = [".+_strelka"]