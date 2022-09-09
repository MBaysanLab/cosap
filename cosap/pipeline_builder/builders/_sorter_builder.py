from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats
from ..._pipeline_config import PipelineKeys, SortingKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Sorter(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    params: Dict = None
    name: str = None
    key: str = PipelineKeys.SORTING

    def __post_init__(self):
        if self.name is None:
            self.name = self.input_step.name

    def _create_config(self) -> Dict:
        filename = self.input_step.get_output()
        output_filename = FileFormats.SORTING_OUTPUT.format(identification=self.name)

        config = {
            self.name: {
                SortingKeys.INPUT: filename,
                SortingKeys.OUTPUT: output_filename,
                SortingKeys.PARAMS: self.params,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.SORTING][self.name][SortingKeys.OUTPUT]

    def get_config(self) -> Dict:
        sorter_config = self._create_config()
        config = {PipelineKeys.SORTING: sorter_config}
        return config
