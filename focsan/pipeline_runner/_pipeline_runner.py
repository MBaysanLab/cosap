from typing import Dict, List

from ..mappers import MapperFactory
from ..variant_callers import VariantCallerFactory
from ._keys import MappingKeys, PipelineKeys


class PipelineRunner:
    def map(self, mapping_configs: List):
        for config in mapping_configs:
            mapper = MapperFactory.create(config[MappingKeys.LIBRARY])
            mapper.map(config)

    def sort(self, sorting_configs: List):
        pass

    def call_variants(self, variant_configs: List):
        for config in variant_configs:
            caller = VariantCallerFactory.create(config[VariantCallingKeys.LIBRARY])
            caller.call_variants(config)

    def run_pipeline(self, pipeline_config: Dict):
        self.validate_pipeline_config(pipeline_config)
        self.map(pipeline_config[PipelineKeys.MAPPING])
        self.sort(pipeline_config[PipelineKeys.SORTING])
        self.call_variants(pipeline_config[PipelineKeys.VARIANT_CALLING])
