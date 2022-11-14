from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import MergingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Merger(_IPipelineStep, _PipelineStep):
    input_step: List[_PipelineStep]
    name: str = None
    key: str = PipelineKeys.MERGE
    next_step: _PipelineStep = None

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

        self.input_step.next_step = self

    def _create_config(self) -> Dict:
        files = []
        for inp in self.input_step:
            files.append(inp.get_output())

        output_filename = FileFormats.MERGING_OUTPUT.format(identification=self.name)

        config = {
            MergingKeys.INPUT: files,
            MergingKeys.OUTPUT: join_paths(
                OutputFolders.PREPROCESSOR, self.key, output_filename
            ),
            MergingKeys.OUTPUT_DIR: join_paths(OutputFolders.PREPROCESSOR, self.key),
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.MERGE][MergingKeys.OUTPUT]

    def get_config(self) -> Dict:
        merger_config = self._create_config()
        config = {PipelineKeys.MERGE: merger_config}
        return config
