from typing import List, Dict
import os
from cosap.mappers import MapperFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import MappingKeys, PipelineKeys


rule mapper:
    input:
        fastqfiles=lambda wildcards: expand(
            config[PipelineKeys.MAPPING][wildcards.identification][
                MappingKeys.INPUT
            ].values()
        ),
    output:
        bam=FileFormats.MAPPING_OUTPUT,
    run:
        mapper = MapperFactory.create(
            config[PipelineKeys.MAPPING][wildcards.identification][MappingKeys.LIBRARY]
        )
        mapper.map(config[PipelineKeys.MAPPING][wildcards.identification])
