from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import CNVCallingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class CNVCaller(_IPipelineStep, _PipelineStep):
    library: str
    normal: _PipelineStep
    tumor: _PipelineStep
    name: str = None
    bed_file: str = None
    key: str = PipelineKeys.CNV
    next_step: _PipelineStep = None

    def __post_init__(self):
        if self.name is None:
            name_temp = []
            if self.normal:
                name_temp.append(self.normal.name)
            if self.tumor:
                name_temp.append(self.tumor.name)
            name_temp.append(self.library)

            self.name = "_".join(name_temp)

        if self.normal:
            self.normal.next_step = self
        if self.tumor:
            self.tumor.next_step = self

    def _create_config(self) -> Dict:
        tumor_input = self.tumor.get_output()
        normal_input = self.normal.get_output()

        output = FileFormats.CNV_OUTPUT.format(identification=self.name)
        output_dir = join_paths(OutputFolders.CNV, self.library)
        config = {
            self.name: {
                CNVCallingKeys.LIBRARY: self.library,
                CNVCallingKeys.NORMAL_INPUT: normal_input,
                CNVCallingKeys.TUMOR_INPUT: tumor_input,
                CNVCallingKeys.OUTPUT_DIR: output_dir,
                CNVCallingKeys.OUTPUT: join_paths(
                    OutputFolders.CNV, self.library, output
                ),
            }
        }
        if self.bed_file:
            config[self.name][CNVCallingKeys.BED_FILE] = self.bed_file
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][CNVCallingKeys.OUTPUT]

    def get_config(self) -> Dict:
        cnv_caller_config = self._create_config()
        config = {self.key: cnv_caller_config}
        return config
