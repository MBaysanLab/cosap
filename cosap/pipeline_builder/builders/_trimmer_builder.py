from dataclasses import dataclass
from typing import Dict, List
from collections import defaultdict

from ..._formats import FileFormats
from ..._pipeline_config import TrimmingKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Trimmer(_IPipelineStep, _PipelineStep):
    reads: List[_PipelineStep]
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self._get_name()

    def _create_config(self) -> Dict:

        read_filenames = {}
        for reader in self.reads:
            read_filenames[reader.read] = reader.get_output()

        if set(read_filenames.keys()) != set(
            map(str, range(1, len(read_filenames) + 1))
        ):
            raise Exception(
                "Inconsistent read numbers, reads should range from 1 to n."
            )

        output_filenames = {}
        for reader in self.reads:
            output_filenames[reader.read] = FileFormats.TRIMMING_OUTPUT.format(
                d=defaultdict(str, identification=self.name, read_no=f"_{reader.read}")
            )

        config = {
            self.name: {
                TrimmingKeys.INPUT: read_filenames,
                TrimmingKeys.OUTPUT: output_filenames,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.TRIM][self.name][TrimmingKeys.OUTPUT]

    def get_config(self) -> Dict:
        trim_config = self._create_config()
        config = {PipelineKeys.TRIM: trim_config}
        return config
