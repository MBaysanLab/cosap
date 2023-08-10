from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import MDUPKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class MDUP(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    name: str = None
    spark: bool = True
    duplicate_handling_method: str = "delete"
    key: str = PipelineKeys.MDUP
    next_step: _PipelineStep = None

    def __post_init__(self):
        self.key = PipelineKeys.MDUP
        if self.name is None:
            self.name = self.input_step.name

        self.input_step.next_step = self

    def _create_config(self) -> Dict:
        output_filename = FileFormats.MDUP_OUTPUT.format(
            identification=self.input_step.name
        )

        config = {
            self.name: {
                MDUPKeys.INPUT: self.input_step.get_output(),
                MDUPKeys.OUTPUT: join_paths(
                    OutputFolders.PREPROCESSOR, self.key, output_filename
                ),
                MDUPKeys.SPARK: self.spark,
                MDUPKeys.DUPLICATE_HANDLING_METHOD: self.duplicate_handling_method,
                MDUPKeys.OUTPUT_DIR: join_paths(OutputFolders.PREPROCESSOR, self.key),
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][MDUPKeys.OUTPUT]

    def get_config(self) -> Dict:
        mdup_config = self._create_config()
        config = {self.key: mdup_config}
        return config
