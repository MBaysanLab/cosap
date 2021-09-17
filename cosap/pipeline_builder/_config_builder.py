from __future__ import annotations
from dataclasses import make_dataclass

from datetime import datetime
from typing import Dict, List
from uuid import uuid4

from .._formats import FileFormats
from .._pipeline_config import (
    BaseRecalibratorKeys,
    IndexingKeys,
    MappingKeys,
    MergingKeys,
    PipelineKeys,
    SortingKeys,
    VariantCallingKeys,
)
from .._version import version
from .builders import _PipelineStep
from .._config import AppConfig


class Pipeline:
    def __init__(self):
        self._pipeline_steps = []

    def _create_config(self):
        config = {
            PipelineKeys.LIBRARY_PATH: AppConfig.LIBRARY_PATH,
            PipelineKeys.CREATION_DATE: datetime.now().strftime(r"%Y-%m-%d %H:%M:%S"),
            PipelineKeys.VERSION: version,
            PipelineKeys.TRIM: dict(),
            PipelineKeys.MAPPING: dict(),
            PipelineKeys.SORTING: dict(),
            PipelineKeys.INDEX: dict(),
            PipelineKeys.MERGE: dict(),
            PipelineKeys.MDUP: dict(),
            PipelineKeys.CALIBRATE: dict(),
            PipelineKeys.VARIANT_CALLING: dict(),
            PipelineKeys.FINAL_OUTPUT: list(),
        }
        return config

    def add(self, step: _PipelineStep) -> Pipeline:
        self._pipeline_steps.append(step)

        return self

    def build(self) -> Dict:
        pipeline_config = self._create_config()

        for step in self._pipeline_steps:
            step_config = step.get_config()
            step_output = step.get_output()

            if type(step_output) == str:
                pipeline_config[PipelineKeys.FINAL_OUTPUT].append(step_output)
            elif type(step_output) == list:
                pipeline_config[PipelineKeys.FINAL_OUTPUT].extend(step_output)
            elif type(step_output) == dict:
                pipeline_config[PipelineKeys.FINAL_OUTPUT].extend(step_output.values())

            for key, values in step_config.items():
                for k, v in values.items():
                    pipeline_config[key][k] = v

        # TODO: insert validation here
        return pipeline_config
