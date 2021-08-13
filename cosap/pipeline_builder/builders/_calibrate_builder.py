from dataclasses import dataclass
from typing import Dict

from .._formats import FileFormats
from .._pipeline_config import BaseRecalibratorKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Recalibrator(_IPipelineStep, _PipelineStep):
    input: str
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self._get_name()

    def _create_config(self) -> Dict:
        filename = self.input.get_output()
        output_filename = FileFormats.MERGING_OUTPUT.format(identification=self.name)
        table_filename = FileFormats.CALIBRATION_TABLE.format(identification=self.name)

        config = {
            BaseRecalibratorKeys.INPUT: filename,
            BaseRecalibratorKeys.TABLE: table_filename,
            BaseRecalibratorKeys.OUTPUT: output_filename,
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.CALIBRATE][BaseRecalibratorKeys.OUTPUT]

    def get_config(self) -> Dict:
        calibration_config = self._create_config()
        config = {PipelineKeys.CALIBRATE: calibration_config}
        return config
