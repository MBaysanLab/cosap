from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, SortingKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Sorter(_IPipelineStep, _PipelineStep):
    input: _PipelineStep
    params: Dict
    name: str = None
    key: str = PipelineKeys.SORTING

    def __post_init__(self):
        if self.name is None:
            self.name = self._get_name()

    def _create_config(self) -> Dict:
        filename = self.input.get_output()
        output_filename = FileFormats.SORTING_OUTPUT.format(identification=self.name)

        config = {
            SortingKeys.INPUT: filename,
            SortingKeys.OUTPUT: join_paths(
                OutputFolders.PREPROCESSOR, self.keys, output_filename
            ),
            SortingKeys.OUTPUT_DIR: join_paths(OutputFolders.PREPROCESSOR, self.keys),
            SortingKeys.PARAMS: self.params,
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.SORTING][SortingKeys.OUTPUT]

    def get_config(self) -> Dict:
        sorter_config = self._create_config()
        config = {PipelineKeys.SORTING: sorter_config}
        return config
