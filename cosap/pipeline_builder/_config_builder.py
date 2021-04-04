from datetime import datetime
from typing import Dict, List
from uuid import uuid4

from ._formats import FileFormats
from ._version import version
from ._pipeline_steps import _PipelineStep
from ._mapper_builder import Mapper
from ._variant_builder import VariantCaller
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
    # ID_LENGTH = 6

    def __init__(self):
        # self._id_generator = lambda: uuid4().hex[: self.ID_LENGTH].upper()
        # self._identities = set()
        self._pipeline_steps = []

    def _create_config(self):
        config = {
            PipelineKeys.CREATION_DATE: datetime.now().strftime(
                r"%Y-%m-%d %H:%M:%S"
            ),
            PipelineKeys.VERSION: version,
            PipelineKeys.MAPPING: list(),
            PipelineKeys.SORTING: list(),
            PipelineKeys.INDEX: list(),
            PipelineKeys.MERGE: list(),
            PipelineKeys.CALIBRATE: list(),
            PipelineKeys.VARIANT_CALLING: list(),
        }
        return config

    def add(self, step: _PipelineStep):
        if isinstance(step, Mapper):
            self._mapper_steps.append(step)
        elif isinstance(step, VariantCaller):
            self._vc_steps.append(step)
        else:
            raise Exception(f"Pipeline step {type(step)} is unknown.")
        
    def build(self) -> Dict:
        pipeline_config = self._create_config()
        connections = {}
        configs = []

        for mapper in self._mapper_steps:
            connections[mapper.output] = mapper.get_output_file_name()
            mapper_config = mapper.get_config()
            configs.append(mapper_config)

        for caller in self._vc_steps:
            caller.set_inputs(
                tumor=connections.get(caller.tumor, caller.tumor_filename),
                germline=connections.get(caller.germline, caller.germline_filename),
            )

            connections[caller.output] = caller.get_output_file_name()
            caller_config = caller.get_config()
            configs.append(caller_config)

        for step_config in configs:
            for step, config in step_config.items():
                pipeline_config[step].extend(config)

        return pipeline_config


# pipeline = Pipeline()

# # read file
# germline_files = [
#     "/mount/data/sample_1/germline_1.fastq",
#     "/mount/data/sample_1/germline_2.fastq",
# ]

# tumor_files = [
#     "/mount/data/sample_1/tumor_1.fastq",
#     "/mount/data/sample_1/tumor_2.fastq",
# ]

# # add mapping step
# pipeline.add(Mapper(
#     library="bwa",
#     read1=germline_files[0],
#     read2=germline_files[1],
#     params=params,
#     output="mapped-germline",
# ))

# pipeline.add(Mapper(
#     library="bwa",
#     read1=tumor_files[0],
#     read2=tumor_files[1],
#     params=params,
#     output="mapped-tumor",
# ))

# # add variant calling
# pipeline.add(VariantCaller(
#     library="varscan", 
#     germline="mapped-germline", 
#     tumor="mapped-tumor", 
#     params=params, 
#     output="vcf-1",
# ))

# # add variant annotation
# pipeline.add(VariantAnnotation(
#     library="annovar",
#     params=params,
#     input="vcf-1",
#     output="annotated-1"
# ))

# pipeline_config = pipeline.build()

