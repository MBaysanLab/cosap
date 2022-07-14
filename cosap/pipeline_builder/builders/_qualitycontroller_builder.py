from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, FolderFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, QualityControlKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class QualityContoller(_IPipelineStep, _PipelineStep):
    input_step: List[_PipelineStep]
    bed_file: str = None
    name: str = None
    key: str = PipelineKeys.QUALITY_CONTROL

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        filename = self.input.get_output()
        raw_output_foldername: FolderFormats.QUALITY_CONTROLLER_OUTPUT.format(
            identification=self.name
        )
        output_filename = FileFormats.QUALIMAP_HTML_OUTPUT.format(
            identification=self.name
        )

        config = {
            self.name: {
                QualityControlKeys.INPUT: filename,
                QualityControlKeys.RAW_OUTPUT: raw_output_foldername,
                QualityControlKeys.OUTPUT: join_paths(
                    OutputFolders.BAMQC, output_filename
                ),
                QualityControlKeys.OUTPUT_DIR: OutputFolders.BAMQC,
            }
        }
        if self.bed_file is not None:
            config[self.name][QualityControlKeys.BED_FILE] = self.bed_file
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.QUALITY_CONTROL][QualityControlKeys.OUTPUT]

    def get_config(self) -> Dict:
        quality_controller_config = self._create_config()
        config = {PipelineKeys.QUALITY_CONTROL: quality_controller_config}
        return config
