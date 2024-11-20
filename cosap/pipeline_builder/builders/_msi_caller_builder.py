from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import MSICallingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class MSICaller(_IPipelineStep, _PipelineStep):
    library: str
    normal: _PipelineStep
    tumor: _PipelineStep
    name: str = None
    key: str = PipelineKeys.MSI
    next_step: _PipelineStep = None
    output_dir: str = None

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

        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.MSI, self.library)

    def _create_config(self) -> Dict:
        tumor_input = self.tumor.get_output()
        normal_input = self.normal.get_output()

        output = FileFormats.MSI_OUTPUT.format(identification=self.name)
        config = {
            self.name: {
                MSICallingKeys.LIBRARY: self.library,
                MSICallingKeys.NORMAL_INPUT: normal_input,
                MSICallingKeys.TUMOR_INPUT: tumor_input,
                MSICallingKeys.OUTPUT: output,
                MSICallingKeys.OUTPUT_DIR: self.output_dir,
                MSICallingKeys.LOG_FILE: self.log_file,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(
            self.output_dir, config[self.key][self.name][MSICallingKeys.OUTPUT]
        )

    def get_config(self) -> Dict:
        msi_caller_config = self._create_config()
        config = {self.key: msi_caller_config}
        return config
