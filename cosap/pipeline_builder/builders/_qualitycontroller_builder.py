from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, QualityControlKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class QualityController(_IPipelineStep, _PipelineStep):
    library: str
    input_step: List[_PipelineStep]
    bed_file: str = None
    name: str = None
    key: str = PipelineKeys.QUALITY_CONTROL
    next_step: _PipelineStep = None
    output_dir: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

        self.input_step.next_step = self

        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.BAMQC, self.library)

    def _create_config(self) -> Dict:
        input_filename = self.input_step.get_output()

        if self.library.lower() == "qualimap":
            raw_output_dir = join_paths(OutputFolders.BAMQC, self.library, self.name)
            config = {
                self.name: {
                    QualityControlKeys.LIBRARY: self.library,
                    QualityControlKeys.INPUT: input_filename,
                    QualityControlKeys.OUTPUT: raw_output_dir,
                }
            }
        elif self.library.lower() == "mosdepth":
            output_filename = FileFormats.MOSDEPTH_OUTPUT.format(
                identification=self.name
            )
            config = {
                self.name: {
                    QualityControlKeys.LIBRARY: self.library,
                    QualityControlKeys.INPUT: input_filename,
                    QualityControlKeys.OUTPUT: output_filename,
                    QualityControlKeys.OUTPUT_DIR: self.output_dir,
                }
            }

        if self.bed_file is not None:
            config[self.name][QualityControlKeys.BED_FILE] = self.bed_file
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(self.output_dir, config[self.name][QualityControlKeys.OUTPUT])

    def get_config(self) -> Dict:
        quality_controller_config = self._create_config()
        config = {self.key: quality_controller_config}
        return config
