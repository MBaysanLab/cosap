from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats
from ..._pipeline_config import IndexingKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Indexer(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    name: str = None
    key: str = PipelineKeys.INDEX

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _get_file_prefix(self, filename):
        return filename.split("_")[0]

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        prefix = self._get_file_prefix(filename)
        output_filename = FileFormats.INDEXING_OUTPUT.format(
            prefix=prefix, identification=self.name
        )

        config = {
            self.name: {
                IndexingKeys.INPUT: filename,
                IndexingKeys.OUTPUT: output_filename,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.INDEX][self.name][IndexingKeys.OUTPUT]

    def get_config(self) -> Dict:
        indexer_config = self._create_config()
        config = {PipelineKeys.INDEX: indexer_config}
        return config
