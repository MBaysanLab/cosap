from __future__ import annotations

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
            PipelineKeys.MAPPING: dict(),
            PipelineKeys.SORTING: dict(),
            PipelineKeys.INDEX: dict(),
            PipelineKeys.MERGE: dict(),
            PipelineKeys.CALIBRATE: dict(),
            PipelineKeys.VARIANT_CALLING: dict(),
        }
        return config

    def add(self, step: _PipelineStep) -> Pipeline:
        self._pipeline_steps.append(step)
        return self

    def build(self) -> Dict:
        pipeline_config = self._create_config()

        for step in self._pipeline_steps:
            step_config = step.get_config()

            for key, values in step_config.items():
                for k, v in values.items():
                    pipeline_config[key][k] = v

        # TODO: insert validation here
        return pipeline_config


# pipeline = Pipeline()

# # read file
# germline_files = [
#     FastqReader(
#         "/mount/data/sample_1/fastq/germline_1.fastq", platform="illumina", read=1
#     ),
#     FastqReader(
#         "/mount/data/sample_1/fastq/germline_2.fastq", platform="illumina", read=2
#     ),
# ]

# tumor_files = [
#     FastqReader(
#         "/mount/data/sample_1/fastq/tumor_1.fastq", platform="illumina", read=1
#     ),
#     FastqReader(
#         "/mount/data/sample_1/fastq/tumor_2.fastq", platform="illumina", read=2
#     ),
# ]

# # add mapping step
# mapper_1 = Mapper(
#     library="bwa",
#     reads=germline_files,
#     params=params,
# )

# mapper_2 = Mapper(
#     library="bwa",
#     reads=tumor_files,
#     params=params,
# )

# # add variant calling
# caller_1 = VariantCaller(
#     library="varscan",
#     germline=mapper_1,
#     tumor=mapper_2,
#     params=params,
# )

# germline_bam = BamReader("/mount/data/sample_1/bam/germline.bam")
# tumor_bam = BamReader("/mount/data/sample_1/bam/tumor.bam")

# caller_2 = VariantCaller(
#     library="varscan",
#     germline=germline_bam,
#     tumor=tumor_bam,
#     params=params,
# )

# pipeline = (
#     Pipeline()
#     .add(mapper_1)
#     .add(mapper_2)
#     .add(caller_1)
#     .add(caller_2)
# )

# pipeline_config = pipeline.build()
