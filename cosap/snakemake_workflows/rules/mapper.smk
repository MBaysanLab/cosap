from typing import List, Dict
import os
from cosap.tools.mappers import MapperFactory
from cosap._formats import FolderedOutputs
from cosap._pipeline_config import MappingKeys, PipelineKeys


rule mapper:
    input:
        fastqfiles=lambda wildcards: expand(
            config[PipelineKeys.MAPPING][wildcards.identification][
                MappingKeys.INPUT
            ].values()
        ),
    output:
        bam=FolderedOutputs.MAPPING_OUTPUT,
    run:
        mapper = MapperFactory.create(
            config[PipelineKeys.MAPPING][wildcards.identification][MappingKeys.LIBRARY]
        )
        mapper.map(config[PipelineKeys.MAPPING][wildcards.identification], config["device"])
