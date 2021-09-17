from dataclasses import dataclass
from os import W_OK


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
    VARIANT_CALLING: str = "variant-calling"
    FINAL_OUTPUT: str = "final-output"


@dataclass
class PipelineBaseKeys:
    LIBRARY: str = "library"
    PARAMS: str = "params"
    INPUT: str = "input"
    OUTPUT: str = "output"


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
    PAIRED_INPUT_1: str = "input_1"
    PAIRED_INPUT_2: str = "input_2"
    PAIRED_OUTPUT_1: str = "output_1"
    PAIRED_OUTPUT_2: str = "output_2"


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
