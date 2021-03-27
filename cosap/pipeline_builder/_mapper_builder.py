from ._formats import FileFormats
from .pipeline_config import (
    BaseRecalibratorKeys,
    IndexingKeys,
    MappingKeys,
    MergingKeys,
    PipelineKeys,
    SortingKeys,
)

from ._pipeline_steps import _PipelineStep


class Mapper(_PipelineStep):
    _GERMLINE = "germline"
    _TUMOR = "tumor"

    def __init__(
        self,
        library: str,
        input_files: List,
        is_tumor: bool,
        library_params: Dict,
        name: str = None,
    ):
        self.name = name
        mapping_config = self._create_mapping_config(
            library=library,
            input_files=input_files,
            library_params=library_params,
        )
        # TODO: indexing is done multiple times, this isnt final
        sorting_config = self._create_sorting_config(mapping_config)
        index_config = self._create_index_config(sorting_config)
        merge_config = self._create_merge_config(index_config)
        calibrate_config = self._create_calibrate_config(merge_config)

        self.steps = {
            PipelineKeys.MAPPING: mapping_config,
            PipelineKeys.SORTING: sorting_config,
            PipelineKeys.INDEX: index_config,
            PipelineKeys.MERGE: merge_config,
            PipelineKeys.CALIBRATE: calibrate_config,
        }

    def _create_mapping_config(
        self,
        library: str,
        input_files: List,
        is_tumor: bool,
        library_params: Dict,
    ) -> Dict:
        output_filename = FileFormats.MAPPING_OUTPUT.format(identification=name)
        # TODO: if there are only two options, convert this to boolean
        sample_type = self._GERMLINE
        if is_tumor:
            sample_type = self._TUMOR

        config = {
            MappingKeys.LIBRARY: library,
            MappingKeys.INPUTS: input_files,
            MappingKeys.OUTPUT: output_filename,
            MappingKeys.SAMPLE_TYPE: sample_type,
            MappingKeys.PARAMS: library_params,
        }
        return config

    def _create_sorting_config(self, mapping_config: Dict) -> Dict:
        output_filename = FileFormats.SORTING_OUTPUT.format(name=name)
        config = {
            SortingKeys.INPUT: mapping_config[MappingKeys.OUTPUT],
            SortingKeys.OUTPUT: output_filename,
            SortingKeys.PARAMS: {},
        }
        return config

    def _create_index_config(self, sorting_config: Dict) -> Dict:
        output_filename = FileFormats.INDEXING_OUTPUT.format(name=name)
        config = {
            IndexingKeys.INPUT: sorting_config[SortingKeys.OUTPUT],
            IndexingKeys.OUTPUT: output_filename,
            IndexingKeys.PARAMS: {},
        }
        return config

    def _create_merge_config(self, index_config: Dict) -> Dict:
        output_filename = FileFormats.MERGING_OUTPUT.format(name=name)
        # TODO: this is probably gonna be multiple files
        input_files = []
        config = {
            MergingKeys.INPUTS: input_files,
            MergingKeys.OUTPUT: output_filename,
        }
        return config

    def _create_calibrate_config(self, merge_config: Dict) -> Dict:
        output_filename = FileFormats.MERGING_OUTPUT.format(name=name)
        table_filename = FileFormats.CALIBRATION_TABLE.format(name=name)
        config = {
            BaseRecalibratorKeys.INPUT: merge_config[MergingKeys.OUTPUT],
            BaseRecalibratorKeys.TABLE: table_filename,
            BaseRecalibratorKeys.OUTPUT: output_filename,
        }
        return config
