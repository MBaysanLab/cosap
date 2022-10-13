from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, TrimmingKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Trimmer(_IPipelineStep, _PipelineStep):
    reads: List[_PipelineStep]
    name: str = None
    key: str = PipelineKeys.TRIM

    def __post_init__(self):
        self.key = PipelineKeys.TRIM
        if self.name == None:
            self.name = "_".join(set(step.name for step in self.reads[::-1]))

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
            output_filename = FileFormats.TRIMMING_OUTPUT.format(
                identification=reader.name, pair=reader.read
            )
            output_filenames[reader.read] = join_paths(
                OutputFolders.TRIMMING, output_filename
            )

        report_filename = FileFormats.TRIMMING_REPORT_OUTPUT.format(
            identification=self.reads[0].name
        )
        config = {
            self.name: {
                TrimmingKeys.INPUT: read_filenames,
                TrimmingKeys.OUTPUT: output_filenames,
                TrimmingKeys.REPORT_OUTPUT: join_paths(
                    OutputFolders.TRIMMING, report_filename
                ),
                TrimmingKeys.OUTPUT_DIR: OutputFolders.TRIMMING,
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
