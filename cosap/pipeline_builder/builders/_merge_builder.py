from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats
from ..._pipeline_config import MergingKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Merger(_IPipelineStep, _PipelineStep):
    input_step: List[_PipelineStep]
    name: str = None
    key: str = PipelineKeys.MERGE

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        files = []
        for inp in self.inputs:
            files.append(inp.get_output())

        output_filename = FileFormats.MERGING_OUTPUT.format(identification=self.name)

        config = {
            MergingKeys.INPUTS: files,
            MergingKeys.OUTPUT: output_filename,
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.MERGE][MergingKeys.OUTPUT]

    def get_config(self) -> Dict:
        merger_config = self._create_config()
        config = {PipelineKeys.MERGE: merger_config}
        return config
