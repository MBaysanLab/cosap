from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import IndexingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Indexer(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    key: str = PipelineKeys.INDEX
    next_step: _PipelineStep = None

    def __post_init__(self):
        self.name = self.input_step.get_output()

        self.input_step.next_step = self

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        output_filename = FileFormats.INDEXING_OUTPUT.format(
            bam_file = filename
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
