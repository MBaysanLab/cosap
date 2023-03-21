from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import BaseRecalibratorKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Recalibrator(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    name: str = None
    key: str = PipelineKeys.CALIBRATE
    next_step: _PipelineStep = None
    bed_file: str = None

    def __post_init__(self):
        self.key = PipelineKeys.CALIBRATE
        if self.name is None:
            self.name = self.input_step.name

        self.input_step.next_step = self

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        output_filename = FileFormats.CALIBRATION_OUTPUT.format(
            identification=self.name
        )
        table_filename = FileFormats.CALIBRATION_TABLE.format(identification=self.name)

        config = {
            self.name: {
                BaseRecalibratorKeys.INPUT: filename,
                BaseRecalibratorKeys.TABLE: join_paths(
                    OutputFolders.CALIBRATION, table_filename
                ),
                BaseRecalibratorKeys.OUTPUT: join_paths(
                    OutputFolders.CALIBRATION, output_filename
                ),
                BaseRecalibratorKeys.OUTPUT_DIR: OutputFolders.CALIBRATION,
            },
        }
        if self.bed_file is not None:
            config[self.name][BaseRecalibratorKeys.BED_FILE] = self.bed_file
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][BaseRecalibratorKeys.OUTPUT]

    def get_config(self) -> Dict:
        calibration_config = self._create_config()
        config = {self.key: calibration_config}

        return config
