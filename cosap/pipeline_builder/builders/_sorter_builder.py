from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, SortingKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Sorter(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    sorting_method: str = "coordinate"
    params: Dict = None
    name: str = None
    key: str = PipelineKeys.SORTING
    next_step: _PipelineStep = None

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

        self.input_step.next_step = self

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        output_filename = FileFormats.SORTING_OUTPUT.format(identification=self.name)

        config = {
            self.name: {
                SortingKeys.INPUT: filename,
                SortingKeys.SORTING_METHOD: self.sorting_method,
                SortingKeys.OUTPUT: output_filename,
                SortingKeys.PARAMS: self.params,
                SortingKeys.LOG_FILE: self.log_file,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][SortingKeys.OUTPUT]

    def get_config(self) -> Dict:
        sorter_config = self._create_config()
        config = {self.key: sorter_config}
        return config
