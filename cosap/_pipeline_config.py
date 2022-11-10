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
    QUALITY_CONTROL: str = "quality_control"
    FINAL_OUTPUT: str = "final-output"
    LOG: str = "log"


@dataclass
class PipelineBaseKeys:
    LIBRARY: str = "library"
    PARAMS: str = "params"
    INPUT: str = "input"
    OUTPUT: str = "output"
    SNAKEMAKE_OUTPUT: str = "snakemake_output"
    OUTPUT_DIR: str = "output_dir"


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
class QualityControlKeys(PipelineBaseKeys):
    BED_FILE: str = "bed"
    RAW_OUTPUT: str = "raw_output"

@dataclass
class SortingKeys(PipelineBaseKeys):
    BAM_DIR: str = "bam-dir"


@dataclass
class IndexingKeys(PipelineBaseKeys):
    BAM_DIR: str = "bam-dir"


@dataclass
class SplitKeys(PipelineBaseKeys):
    pass


@dataclass
class MergingKeys(PipelineBaseKeys):
    pass


@dataclass
class MDUPKeys(PipelineBaseKeys):
    METRICS: str = "metrics"


@dataclass
class BaseRecalibratorKeys(PipelineBaseKeys):
    TABLE: str = "table"


@dataclass
class VariantCallingKeys(PipelineBaseKeys):
    GERMLINE_INPUT: str = "normal_input"
    TUMOR_INPUT: str = "tumor_input"
    GERMLINE_SAMPLE_NAME: str = "germline_sample_name"
    TUMOR_SAMPLE_NAME: str = "tumor_sample_name"
    UNFILTERED_VARIANTS_OUTPUT: str = "unfiltered_variants"
    FILTERED_VARIANTS_OUTPUT: str = "filtered_variants"
    GVCF_OUTPUT: str = "gvcf_output"
    SNP_OUTPUT: str = "snp_output"
    INDEL_OUTPUT: str = "indel_output"
    OTHER_VARIANTS_OUTPUT: str = "other_variants_output"
    PILEUPS: str = "pileups"


@dataclass
class AnnotatorKeys(PipelineBaseKeys):
    AVOUTPUT: str = "av_output"


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
class SnakemakeConstraints:
    PY2_VARIANT_CALLERS: Pattern = ".+_strelka"
