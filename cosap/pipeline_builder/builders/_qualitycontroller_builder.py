from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, FolderFormats
from ..._pipeline_config import PipelineKeys, QualityControlKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from ..._utils import join_paths


@dataclass
class QualityContoller(_IPipelineStep, _PipelineStep):
    library: str 
    input_step: List[_PipelineStep]
    bed_file: str = None
    name: str = None
    key: str = PipelineKeys.QUALITY_CONTROL

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        
        if self.library.lower() == "qualimap":
            raw_output_folderdir = FolderFormats.QUALIMAP_OUTPUT.format(
                identification=self.name
            )
            output_filename = FileFormats.QUALIMAP_PDF_OUTPUT.format(
                identification=self.name
            )
            config = {
            self.name: {
                QualityControlKeys.LIBRARY: self.library,
                QualityControlKeys.INPUT: filename,
                QualityControlKeys.RAW_OUTPUT: raw_output_folderdir,
                QualityControlKeys.OUTPUT: join_paths(raw_output_folderdir, output_filename),
            }
        }
        elif self.library.lower() == "mosdepth":
            output_filename = FileFormats.MOSDEPTH_OUTPUT.format(
                identification=self.name
            )
            config = {
            self.name: {
                QualityControlKeys.LIBRARY: self.library,
                QualityControlKeys.INPUT: filename,
                QualityControlKeys.OUTPUT: output_filename,
            }
        }
            
       
        if self.bed_file is not None:
            config[self.name][QualityControlKeys.BED_FILE] = self.bed_file
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][QualityControlKeys.OUTPUT]

    def get_config(self) -> Dict:
        quality_controller_config = self._create_config()
        config = {PipelineKeys.QUALITY_CONTROL: quality_controller_config}
        return config
