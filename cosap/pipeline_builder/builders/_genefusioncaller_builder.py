from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import GeneFusionCallingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class GeneFusionCaller(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    library: str
    name: str = None
    key: str = PipelineKeys.GENEFUSION
    next_step: _PipelineStep = None
    output_dir: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = "_".join(set(step.name for step in self.input_step[::-1]))

        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.GENE_FUSION, self.library)

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

        json_output = FileFormats.GENEFUSION_OUTPUT.format(identification=self.name)

        config = {
            self.name: {
                GeneFusionCallingKeys.LIBRARY: self.library,
                GeneFusionCallingKeys.INPUT: read_filenames,
                GeneFusionCallingKeys.OUTPUT: json_output,
                GeneFusionCallingKeys.OUTPUT_DIR: self.output_dir,
                GeneFusionCallingKeys.LOG_FILE: self.log_file,
            }
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(
            self.output_dir, config[self.key][self.name][GeneFusionCallingKeys.OUTPUT]
        )

    def get_config(self) -> Dict:
        genefusion_caller_config = self._create_config()
        config = {self.key: genefusion_caller_config}
        return config
