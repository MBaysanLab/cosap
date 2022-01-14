from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats
from ..._pipeline_config import ElprepKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Elprep(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    name: str = None

    def __post_init__(self):
        self.key = PipelineKeys.ELPREP_PROCESS
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        output_filename = FileFormats.ELPREP_CALIBRATION_OUTPUT.format(
            identification=self.name
        )
        table_filename = FileFormats.CALIBRATION_TABLE.format(identification=self.name)
        config = {
            self.name: {
                ElprepKeys.INPUT: self.input_step.get_output(),
                ElprepKeys.TABLE: table_filename,
                ElprepKeys.OUTPUT: output_filename,
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][ElprepKeys.OUTPUT]

    def get_config(self) -> Dict:
        elprep_config = self._create_config()
        config = {self.key: elprep_config}
        return config
