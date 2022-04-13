from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats
from ..._pipeline_config import BaseRecalibratorKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Recalibrator(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    name: str = None
    key: str = PipelineKeys.CALIBRATE

    def __post_init__(self):
        self.key = PipelineKeys.CALIBRATE
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        output_filename = FileFormats.CALIBRATION_OUTPUT.format(
            identification=self.name
        )
        table_filename = FileFormats.CALIBRATION_TABLE.format(identification=self.name)

        config = {
            self.name: {
                BaseRecalibratorKeys.INPUT: filename,
                BaseRecalibratorKeys.TABLE: table_filename,
                BaseRecalibratorKeys.OUTPUT: output_filename,
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][BaseRecalibratorKeys.OUTPUT]

    def get_config(self) -> Dict:
        calibration_config = self._create_config()
        config = {self.key: calibration_config}
        return config
