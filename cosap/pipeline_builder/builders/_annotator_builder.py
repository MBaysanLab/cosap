from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats
from ..._pipeline_config import MDUPKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Annotator(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        output_filename = FileFormats.MDUP_OUTPUT.format(
            identification=self.input_step.name
        )

        config = {
            self.name: {
                MDUPKeys.INPUT: self.input_step.get_output(),
                MDUPKeys.OUTPUT: output_filename,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.MDUP][self.name][MDUPKeys.OUTPUT]

    def get_config(self) -> Dict:
        mdup_config = self._create_config()
        config = {PipelineKeys.MDUP: mdup_config}
        return config
