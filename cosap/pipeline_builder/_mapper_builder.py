from ._formats import FileFormats
from typing import List, Dict
from .pipeline_config import (
    BaseRecalibratorKeys,
    IndexingKeys,
    MappingKeys,
    MergingKeys,
    PipelineKeys,
    SortingKeys,
)

from ._pipeline_steps import _PipelineStep
from dataclasses import dataclass
from copy import copy


@dataclass
class Mapper(_PipelineStep):
    library: str
    read1: str
    read2: str
    params: Dict
    output: str

    def _create_mapping_config(self) -> Dict:
        output_filename = FileFormats.MAPPING_OUTPUT.format(identification=self.name)

        config = {
            MappingKeys.LIBRARY: self.library,
            MappingKeys.READ1: self.read1,
            MappingKeys.READ2: self.read2
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

    def get_output_file_name(self):
        return self.config[PipelineKeys.CALIBRATE][BaseRecalibratorKeys.OUTPUT]

    def get_config(self):
        return copy(self.config)

    @property
    def config(self):
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
    
