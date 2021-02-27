from .._pipeline_config import PipelineConfig
from .sorter import SorterFactory


class PreProcessing:
    @classmethod
    def preprocess(cls, pipeline_config: PipelineConfig):
        BAMSorter.sort(pipeline_config=pipeline_config)
