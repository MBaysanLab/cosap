from dataclasses import dataclass
from typing import Pattern


@dataclass
class PipelineKeys:
    WORKDIR: str = "workdir"
    LIBRARY_PATH: str = "library_path"
    CREATION_DATE: str = "creation-date"
    VERSION: str = "cosap_version"
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
    GENEFUSION: str = "genefusion"
    MSI: str = "msi"
    CNV: str = "cnv"


@dataclass
class PipelineBaseKeys:
    LIBRARY: str = "library"
    PARAMS: str = "params"
    INPUT: str = "input"
    OUTPUT: str = "output"
    SNAKEMAKE_OUTPUT: str = "snakemake_output"
    OUTPUT_DIR: str = "output_dir"
    BED_FILE: str = "bed"


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
    REPORT_OUTPUT: str = "report"


@dataclass
class QualityControlKeys(PipelineBaseKeys):
    RAW_OUTPUT: str = "raw_output"


@dataclass
class SortingKeys(PipelineBaseKeys):
    BAM_DIR: str = "bam-dir"
    SORTING_METHOD: str = "sort-type"


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
    SPARK: str = "spark"
    DUPLICATE_HANDLING_METHOD: str = "duplicate_handling_method"


@dataclass
class BaseRecalibratorKeys(PipelineBaseKeys):
    TABLE: str = "table"


@dataclass
class VariantCallingKeys(PipelineBaseKeys):
    OUTPUT_TYPE: str = "output_type"
    GERMLINE_INPUT: str = "normal_input"
    TUMOR_INPUT: str = "tumor_input"
    GERMLINE_SAMPLE_NAME: str = "germline_sample_name"
    TUMOR_SAMPLE_NAME: str = "tumor_sample_name"
    UNFILTERED_VARIANTS_OUTPUT: str = "unfiltered_variants"
    ALL_VARIANTS_OUTPUT: str = "all_variants"
    GVCF_OUTPUT: str = "gvcf_output"
    SNP_OUTPUT: str = "snp_output"
    INDEL_OUTPUT: str = "indel_output"
    OTHER_VARIANTS_OUTPUT: str = "other_variants_output"
    PILEUPS: str = "pileups"


@dataclass
class GeneFusionCallingKeys(PipelineBaseKeys):
    pass


@dataclass
class MSICallingKeys(PipelineBaseKeys):
    NORMAL_INPUT: str = "normal_input"
    TUMOR_INPUT: str = "tumor_input"


@dataclass
class CNVCallingKeys(PipelineBaseKeys):
    NORMAL_INPUT: str = "normal_input"
    TUMOR_INPUT: str = "tumor_input"
    OUTPUT_DIR: str = "output_dir"


@dataclass
class AnnotatorKeys(PipelineBaseKeys):
    AVOUTPUT: str = "av_output"
    INPUT_TYPE: str = "input_type"


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
    PY2_VARIANT_CALLERS: Pattern = ".+_strelka|.+_manta"
