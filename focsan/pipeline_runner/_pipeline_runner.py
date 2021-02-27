from typing import Dict, List
from ._keys import PipelineKeys, MappingKeys
from ..mappers import MapperFactory


class PipelineRunner:
    def map(self, mapping_configs: List):
        for config in mapping_configs:
            mapper = MapperFactory.create(config[MappingKeys.LIBRARY])
            mapper.map(config)

    def sort(self, sorting_configs: List):
        pass

    def call_variants(self, variant_configs: List):
        pass

    def run_pipeline(self, pipeline_config: Dict):
        self.validate_pipeline_config(pipeline_config)
        self.map(pipeline_config[PipelineKeys.MAPPING])
        self.sort(pipeline_config[PipelineKeys.SORTING])
        self.call_variants(pipeline_config[PipelineKeys.VARIANT_CALLING])
