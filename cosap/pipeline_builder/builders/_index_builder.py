from dataclasses import dataclass
from typing import Dict

from .._formats import FileFormats
from .._pipeline_config import IndexingKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Indexer(_IPipelineStep, _PipelineStep):
    input: str
    params: Dict
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self._get_name()

    def _create_config(self) -> Dict:
        filename = self.input.get_output()
        output_filename = FileFormats.INDEXING_OUTPUT.format(identification=self.name)

        config = {
            IndexingKeys.INPUT: filename,
            IndexingKeys.OUTPUT: output_filename,
            IndexingKeys.PARAMS: self.params,
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.INDEX][IndexingKeys.OUTPUT]

    def get_config(self) -> Dict:
        indexer_config = self._create_config()
        config = {PipelineKeys.INDEX: indexer_config}
        return config
