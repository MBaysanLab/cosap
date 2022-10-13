from __future__ import annotations

from dataclasses import make_dataclass
from datetime import datetime
from typing import Dict, List

from .._config import AppConfig
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
            PipelineKeys.ELPREP_PROCESS: dict(),
            PipelineKeys.VARIANT_CALLING: dict(),
            PipelineKeys.ANNOTATION: dict(),
            PipelineKeys.QUALITY_CONTROL: dict(),
            PipelineKeys.FINAL_OUTPUT: list(),
            PipelineKeys.WORKDIR: str,
        }
        return config

    def add(self, step: _PipelineStep) -> Pipeline:
        self._pipeline_steps.append(step)

        return self

    def build(self, workdir: str = None) -> Dict:
        pipeline_config = self._create_config()

        if workdir:
            pipeline_config[PipelineKeys.WORKDIR] = workdir

        # Final output is the outputs of latest added step
        step_key = self._pipeline_steps[-1].key
        for step in self._pipeline_steps:
            step_config = step.get_config()

            for key, values in step_config.items():
                if key == step_key:
                    pipeline_config[PipelineKeys.FINAL_OUTPUT].append(step.get_output())
                for k, v in values.items():
                    pipeline_config[key][k] = v

        # TODO: insert validation here
        return pipeline_config
