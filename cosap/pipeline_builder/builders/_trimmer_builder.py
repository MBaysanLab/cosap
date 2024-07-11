import os
from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import PipelineKeys, TrimmingKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class Trimmer(_IPipelineStep, _PipelineStep):
    input_step: List[_PipelineStep]
    name: str = None
    key: str = PipelineKeys.TRIM
    next_step: _PipelineStep = None
    output_dir = OutputFolders.TRIMMING

    def __post_init__(self):
        self.key = PipelineKeys.TRIM
        if self.name == None:
            self.name = self._create_name_from_input()

    def _create_config(self) -> Dict:
        read_filenames = {}
        for reader in self.input_step:
            read_filenames[reader.read] = reader.get_output()

        if set(read_filenames.keys()) != set(
            map(str, range(1, len(read_filenames) + 1))
        ):
            raise Exception(
                "Inconsistent read numbers, reads should range from 1 to n."
            )

        output_filenames = {}
        for reader in self.input_step:
            output_filename = FileFormats.TRIMMING_OUTPUT.format(
                identification=self.name, pair=reader.read
            )
            output_filenames[reader.read] = output_filename

        report_filename = FileFormats.TRIMMING_REPORT_OUTPUT.format(
            identification=self.input_step[0].name
        )
        config = {
            self.name: {
                TrimmingKeys.INPUT: read_filenames,
                TrimmingKeys.OUTPUT: output_filenames,
                TrimmingKeys.REPORT_OUTPUT: report_filename,
                TrimmingKeys.OUTPUT_DIR: self.output_dir,
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(
            self.output_dir, config[self.key][self.name][TrimmingKeys.OUTPUT]
        )

    def get_config(self) -> Dict:
        trim_config = self._create_config()
        config = {self.key: trim_config}
        return config

    def _create_name_from_input(self):
        """
        Create name from input steps if name is not provided.

        Returns:
            str: Name of the step
        """

        # Get common prefix of the input fastq files
        prefix = os.path.commonprefix(
            [reader.get_output() for reader in self.input_step]
        )

        # Remove trailing underscore and the rest of the name
        prefix = os.path.basename(prefix).rstrip("_")
        return prefix
