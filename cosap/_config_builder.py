from datetime import datetime
from typing import Dict, List
from uuid import uuid4

from ._formats import FileFormats
from ._version import version
from .pipeline_config import MappingKeys, PipelineKeys


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
        self._id_generator = lambda: uuid4().hex[:self.ID_LENGTH].upper()
        self._identities = set()

    def _get_new_id(self) -> str:
        new_id = self._id_generator()
        while new_id in self._identities:
            new_id = self._id_generator()

        self._identities.add(new_id)
        return new_id

    def _create_mapping_config(self, , library: str, input_files: List, is_tumor: bool, library_params: Dict, output_filename: str = None) -> Dict:
        if output_filename is None:
            identification = self._get_new_id()
            output_filename = FileFormats.MAPPING_OUTPUT.format(identification=identification)        

        # TODO: if there are only two options, convert this to boolean
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

    def _create_preprocessing_configs(self, mapping_config: Dict) -> :
        pass

    def add_mapping_step(self, library: str, input_files: List, is_tumor: bool, library_params: Dict, output_filename: str = None):
        mapping_config = self._create_mapping_config(library=library, input_files=input_files, library_params=library_params, output_filename=output_filename)
        # TODO: this could be multiple steps to go along with the other create methods
        self._create_preprocessing_configs(mapping_config=mapping_config)

        self._config[PipelineKeys.MAPPING].append(mapping_config)
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

