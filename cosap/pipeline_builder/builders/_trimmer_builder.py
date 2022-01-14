from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats
from ..._pipeline_config import PipelineKeys, TrimmingKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Trimmer(_IPipelineStep, _PipelineStep):
    reads: List[_PipelineStep]
    name: str = None

    def __post_init__(self):
        self.key = PipelineKeys.TRIM
        if self.name == None:
            self.name = "_".join([step.name for step in self.reads])

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
                identification=reader.name, pair=reader.read
            )

        config = {
            self.name: {
                TrimmingKeys.INPUT: read_filenames,
                TrimmingKeys.OUTPUT: output_filenames,
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][TrimmingKeys.OUTPUT]

    def get_config(self) -> Dict:
        trim_config = self._create_config()
        config = {self.key: trim_config}
        return config
