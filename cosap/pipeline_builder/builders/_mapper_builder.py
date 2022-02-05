from copy import copy
from dataclasses import dataclass
from typing import Dict, List, Union

from ..._formats import FileFormats
from ..._pipeline_config import (
    BaseRecalibratorKeys,
    IndexingKeys,
    MappingKeys,
    MergingKeys,
    PipelineKeys,
    SortingKeys,
)
from ._file_readers import FastqReader
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from ._trimmer_builder import Trimmer


@dataclass
class Mapper(_IPipelineStep, _PipelineStep):
    library: str
    reads: Union[_PipelineStep, List[_PipelineStep]]
    params: Dict
    name: str = None

    def __post_init__(self):
        self.key = PipelineKeys.MAPPING
        if self.name is None:
            if isinstance(self.reads, Trimmer):
                self.name = f"{self.reads.name}_{self.library}"
            else:
                self.name = (
                    "%s" % "-".join(read.name for read in self.reads) + self.library
                )

    def _create_config(self) -> Dict:
        output_filename = FileFormats.MAPPING_OUTPUT.format(identification=self.name)

        if type(self.reads) == Trimmer:
            read_filenames = self.reads.get_output()

        else:
            read_filenames = {}
            for reader in self.reads:
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
                MappingKeys.PARAMS: self.params,
            },
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[self.key][self.name][MappingKeys.OUTPUT]

    def get_config(self) -> Dict:
        mapping_config = self._create_config()
        config = {
            self.key: mapping_config,
        }
        return config
