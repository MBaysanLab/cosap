from dataclasses import dataclass


@dataclass
class PipelineKeys:
    CREATION_DATE: str = "creation-date"
    VERSION: str = "version"
    MAPPING: str = "mapping"
    SORTING: str = "sorting"
    INDEX: str = "index"
    MERGE: str = "merge"
    CALIBRATE: str = "calibrate"
    VARIANT_CALLING: str = "variant-calling"


@dataclass
class PipelineBaseKeys:
    LIBRARY: str = "library"
    PARAMS: str = "params"
    INPUT: str = "input"
    OUTPUT: str = "output"


@dataclass
class MappingKeys(PipelineBaseKeys):
    SAMPLE_TYPE: str = "sample-type"
    READ1: str = "read-1"
    READ2: str = "read-2"
    RG_ID: str = "rg-id"
    RG_SM: str = "rg-sm"
    RG_LB: str = "rg-lb"
    RG_PL: str = "rg-pl"
    RG_PU: str = "rg-pu"


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
    NORMAL_SRC: str = "normal-path"
    TUMOR_SRC: str = "tumor-path"
    OUTPUT_DIR: str = "output-dir"
    PILEUPS: str = "pileups"
