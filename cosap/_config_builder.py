from datetime import datetime
from typing import Dict, List
from uuid import uuid4

from ._formats import FileFormats
from ._version import version
from .pipeline_config import (BaseRecalibratorKeys, IndexingKeys, MappingKeys,
                              MergingKeys, PipelineKeys, SortingKeys)


class ConfigBuilder:
    # TODO: probably just having map and variant call adders
    # should be enough. Preprocessing is always done anyway
    ID_LENGTH = 6

    def __init__(self):
        self._config = {
            PipelineKeys.VERSION: version,
            PipelineKeys.MAPPING: list(),
            PipelineKeys.SORTING: list(),
            PipelineKeys.INDEX: list(),
            PipelineKeys.MERGE: list(),
            PipelineKeys.CALIBRATE: list(),
            PipelineKeys.VARIANT_CALLING: list(),
        }
        self._id_generator = lambda: uuid4().hex[: self.ID_LENGTH].upper()
        self._identities = set()

    def _get_new_id(self) -> str:
        new_id = self._id_generator()
        while new_id in self._identities:
            new_id = self._id_generator()

        self._identities.add(new_id)
        return new_id

    def _create_mapping_config(
        self,
        library: str,
        input_files: List,
        is_tumor: bool,
        library_params: Dict,
        identification: str,
    ) -> Dict:
        # TODO: if there are only two options, convert this to boolean
        output_filename = FileFormats.MAPPING_OUTPUT.format(
            identification=identification
        )
        sample_type = "germline"
        if is_tumor:
            sample_type = "tumor"

        config = {
            MappingKeys.LIBRARY: library,
            MappingKeys.INPUTS: input_files,
            MappingKeys.OUTPUT: output_filename,
            MappingKeys.SAMPLE_TYPE: sample_type,
            MappingKeys.PARAMS: library_params,
        }
        return config

    def _create_sorting_config(self, mapping_config: Dict, identification: str) -> Dict:
        output_filename = FileFormats.SORTING_OUTPUT.format(
            identification=identification
        )
        config = {
            SortingKeys.INPUT: mapping_config[MappingKeys.OUTPUT],
            SortingKeys.OUTPUT: output_filename,
            SortingKeys.PARAMS: {},
        }
        return config

    def _create_index_config(self, sorting_config: Dict, identification: str) -> Dict:
        output_filename = FileFormats.INDEXING_OUTPUT.format(
            identification=identification
        )
        config = {
            IndexingKeys.INPUT: sorting_config[SortingKeys.OUTPUT],
            IndexingKeys.OUTPUT: output_filename,
            IndexingKeys.PARAMS: {},
        }
        return config

    def _create_merge_config(self, index_config: Dict, identification: str) -> Dict:
        output_filename = FileFormats.MERGING_OUTPUT.format(
            identification=identification
        )
        # TODO: this is probably gonna be multiple files
        input_files = []
        config = {
            MergingKeys.INPUTS: input_files,
            MergingKeys.OUTPUT: output_filename,
        }
        return config

    def _create_calibrate_config(self, merge_config: Dict, identification: str) -> Dict:
        output_filename = FileFormats.MERGING_OUTPUT.format(
            identification=identification
        )
        table_filename = FileFormats.CALIBRATION_TABLE.format(
            identification=identification
        )
        config = {
            BaseRecalibratorKeys.INPUT: merge_config[MergingKeys.OUTPUT],
            BaseRecalibratorKeys.TABLE: table_filename,
            BaseRecalibratorKeys.OUTPUT: output_filename,
        }
        return config

    def add_mapping_step(
        self,
        library: str,
        input_files: List,
        is_tumor: bool,
        library_params: Dict,
        identification: str = None,
    ):
        if identification is None:
            # TODO: this id should most likely be regenratable from
            # the inputs, like hash or just str concatenation
            identification = self._get_new_id()

        mapping_config = self._create_mapping_config(
            library=library,
            input_files=input_files,
            library_params=library_params,
            identification=identification,
        )
        sorting_config = self._create_sorting_config(mapping_config, identification)
        index_config = self._create_index_config(sorting_config, identification)
        merge_config = self._create_merge_config(index_config, identification)
        calibrate_config = self._create_calibrate_config(merge_config, identification)

        self._config[PipelineKeys.MAPPING].append(mapping_config)
        self._config[PipelineKeys.SORTING].append(sorting_config)
        self._config[PipelineKeys.INDEX].append(index_config)
        self._config[PipelineKeys.MERGE].append(merge_config)
        self._config[PipelineKeys.CALIBRATE].append(calibrate_config)

        return self

    def add_variant_calling_step(self, config: Dict):
        self._config[PipelineKeys.VARIANT_CALLING].append(config)
        return self

    def build_config(self):
        self._config[PipelineKeys.CREATION_DATE] = datetime.now().strftime(
            r"%Y-%m-%d %H:%M:%S"
        )
        # TODO: should probably do a check here to see if config makes sense
        return self._config
