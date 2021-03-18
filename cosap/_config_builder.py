from datetime import datetime

from ._version import version
from .pipeline_config import PipelineKeys


class ConfigBuilder:
    # TODO: probably just having map and variant call adders
    # should be enough. Preprocessing is always done anyway
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

    def add_map_config(self, config: Dict):
        self._config[PipelineKeys.MAPPING].append(config)
        return self

    def add_sort_config(self, config: Dict):
        self._config[PipelineKeys.SORTING].append(config)
        return self

    def add_index_config(self, config: Dict):
        self._config[PipelineKeys.INDEX].append(config)
        return self

    def add_merge_config(self, config: Dict):
        self._config[PipelineKeys.MERGE].append(config)
        return self

    def add_calibrate_config(self, config: Dict):
        self._config[PipelineKeys.CALIBRATE].append(config)
        return self

    def add_variant_caller_config(self, config: Dict):
        self._config[PipelineKeys.VARIANT_CALLING].append(config)
        return self

    def build_config(self):
        self._config[PipelineKeys.CREATION_DATE] = datetime.now().strftime(
            r"%Y-%m-%d %H:%M:%S"
        )
        # TODO: should probably do a check here to see if config makes sense
        return self._config

