from copy import copy
from dataclasses import dataclass
from typing import Dict, List

from ._formats import FileFormats
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from .pipeline_config import (BaseRecalibratorKeys, IndexingKeys, MappingKeys,
                              MergingKeys, PipelineKeys, SortingKeys)


@dataclass
class Mapper(_IPipelineStep, _PipelineStep):
    library: str
    reads: List
    params: Dict
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self._get_name()

    def _create_mapping_config(self) -> Dict:
        output_filename = FileFormats.MAPPING_OUTPUT.format(identification=self.name)

        read_filenames = {}
        for reader in self.reads:
            read_filenames[reader.read] = reader.get_output()

        if set(read_filenames.keys()) != set(range(1, len(read_filenames))):
            raise Exception(
                "Inconsistent read numbers, reads should range from 1 to n."
            )

        config = {
            MappingKeys.LIBRARY: self.library,
            MappingKeys.INPUT: read_filenames,
            MappingKeys.OUTPUT: output_filename,
            MappingKeys.PARAMS: self.params,
        }
        return config

    def _create_sorting_config(self, mapping_config: Dict) -> Dict:
        output_filename = FileFormats.SORTING_OUTPUT.format(identification=self.name)
        config = {
            SortingKeys.INPUT: mapping_config[MappingKeys.OUTPUT],
            SortingKeys.OUTPUT: output_filename,
            SortingKeys.PARAMS: {},
        }
        return config

    def _create_index_config(self, config: Dict) -> Dict:
        output_filename = FileFormats.INDEXING_OUTPUT.format(identification=self.name)
        config = {
            IndexingKeys.INPUT: config[SortingKeys.OUTPUT],
            IndexingKeys.OUTPUT: output_filename,
            IndexingKeys.PARAMS: {},
        }
        return config

    def _create_merge_config(self, index_config: Dict) -> Dict:
        output_filename = FileFormats.MERGING_OUTPUT.format(identification=self.name)
        # TODO: this is probably gonna be multiple files
        input_files = []
        config = {
            MergingKeys.INPUTS: input_files,
            MergingKeys.OUTPUT: output_filename,
        }
        return config

    def _create_calibrate_config(self, merge_config: Dict) -> Dict:
        output_filename = FileFormats.MERGING_OUTPUT.format(identification=self.name)
        table_filename = FileFormats.CALIBRATION_TABLE.format(identification=self.name)
        config = {
            BaseRecalibratorKeys.INPUT: merge_config[MergingKeys.OUTPUT],
            BaseRecalibratorKeys.TABLE: table_filename,
            BaseRecalibratorKeys.OUTPUT: output_filename,
        }
        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.CALIBRATE][BaseRecalibratorKeys.OUTPUT]

    def get_config(self) -> Dict:
        mapping_config = self._create_mapping_config()
        # TODO: indexing is done multiple times, this isnt final
        sorting_config = self._create_sorting_config(mapping_config)
        index_config = self._create_index_config(sorting_config)
        merge_config = self._create_merge_config(index_config)
        calibrate_config = self._create_calibrate_config(merge_config)

        config = {
            PipelineKeys.MAPPING: mapping_config,
            PipelineKeys.SORTING: sorting_config,
            PipelineKeys.INDEX: index_config,
            PipelineKeys.MERGE: merge_config,
            PipelineKeys.CALIBRATE: calibrate_config,
        }
        return config
