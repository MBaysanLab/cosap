from dataclasses import dataclass, field
from typing import Dict, List, Union

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import MappingKeys, PipelineKeys
from ..._utils import join_paths
from ._file_readers import FastqReader
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from ._trimmer_builder import Trimmer


@dataclass
class Mapper(_IPipelineStep, _PipelineStep):
    library: str
    input_step: Union[_PipelineStep, List[_PipelineStep]]
    params: Dict = field(default_factory=dict)
    name: str = None
    key: str = PipelineKeys.MAPPING
    next_step: _PipelineStep = None
    post_processing: bool = False
    output_dir: str = None

    def __post_init__(self):
        self.key = PipelineKeys.MAPPING
        if self.name is None:
            if isinstance(self.input_step, Trimmer):
                self.name = f"{self.input_step.name}_{self.library}"
            else:
                self.name = (
                    "%s" % "-".join(read.name for read in self.input_step)
                    + self.library
                )

        self.library = self.library.lower()

        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.MAPPING, self.library)

        if isinstance(self.input_step, list):
            for step in self.input_step:
                step.next_step = self
        else:
            self.input_step.next_step = self

        if (MappingKeys.READ_GROUP not in self.params.keys()) or (
            MappingKeys.RG_SM not in self.params[MappingKeys.READ_GROUP].keys()
        ):
            raise Exception("Please specify a sample name for the read group.")

    def _create_config(self) -> Dict:
        output_filename = FileFormats.MAPPING_OUTPUT.format(identification=self.name)

        if type(self.input_step) == Trimmer:
            read_filenames = self.input_step.get_output()

        else:
            read_filenames = {}
            for reader in self.input_step:
                read_filenames[reader.read] = reader.get_output()

        if set(read_filenames.keys()) != set(
            map(str, range(1, len(read_filenames) + 1))
        ):
            raise Exception(
                "Inconsistent read numbers, reads should range from 1 to n."
            )

        config = {
            self.name: {
                MappingKeys.LIBRARY: self.library,
                MappingKeys.INPUT: read_filenames,
                MappingKeys.OUTPUT: output_filename,
                MappingKeys.OUTPUT_DIR: self.output_dir,
                MappingKeys.PARAMS: self.params,
                MappingKeys.POST_PROCESSING: self.post_processing,
                MappingKeys.LOG_FILE: self.log_file,
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(
            self.output_dir, config[self.key][self.name][MappingKeys.OUTPUT]
        )

    def get_config(self) -> Dict:
        mapping_config = self._create_config()
        config = {
            self.key: mapping_config,
        }
        return config
