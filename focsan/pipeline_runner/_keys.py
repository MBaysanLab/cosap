from dataclasses import dataclass

@dataclass
class PipelineKeys:
	CREATION_DATE: str = "creation-date"
	VERSION: str = "version"
	MAPPING: str = "mapping"
	SORTING: str = "sorting"
	VARIANT_CALLING: str = "variant-calling"


@dataclass
class PipelineBaseKeys:
	LIBRARY: str = "library"
	PARAMS: str = "params"


@dataclass
class MappingKeys(PipelineBaseKeys):
	SAMPLE_TYPE: str = "sample-type"


@dataclass
class SortingKeys(PipelineBaseKeys):
	pass


@dataclass
class VariantCallingKeys(PipelineBaseKeys):
	NORMAL_SRC: str = "normal-path"
	TUMOR_SRC: str = "tumor-path"
