from datetime import datetime
from typing import Dict, List
from uuid import uuid4

from ._formats import FileFormats
from ._version import version
from ._pipeline_steps import _PipelineStep
from .pipeline_config import (
    BaseRecalibratorKeys,
    IndexingKeys,
    MappingKeys,
    MergingKeys,
    PipelineKeys,
    SortingKeys,
    VariantCallingKeys,
)


class Pipeline:
    # TODO: probably just having map and variant call adders
    # should be enough. Preprocessing is always done anyway
    # ID_LENGTH = 6

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
        # self._id_generator = lambda: uuid4().hex[: self.ID_LENGTH].upper()
        # self._identities = set()
        self._steps = []

    def add(self, step: _PipelineStep):
        self._steps.append(step)

    def build(self) -> Dict:
        self._config[PipelineKeys.CREATION_DATE] = datetime.now().strftime(
            r"%Y-%m-%d %H:%M:%S"
        )
        # TODO: do a check here for validating each step
        for pipeline_step in self._steps:
            for key, params in pipeline_step.steps.items():
                if key not in self._config:
                    raise Exception(f"{key} is not a known pipeline step.")
                self._config[key].append(params)
        return self._config
