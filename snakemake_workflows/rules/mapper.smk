from typing import List, Dict
import os
from cosap.mappers import MapperFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import MappingKeys,PipelineKeys

rule bwa:
    output:
        bam=FileFormats.MAPPING_OUTPUT.format(id="{identification}", mapper="bwa")
    run:
        bwa_mapper = MapperFactory.create(config[PipelineKeys.MAPPING][wildcards.identification][MappingKeys.LIBRARY])
        bwa_mapper.map(config[PipelineKeys.MAPPING][wildcards.identification])

rule bowtie:
    output:
        bam=FileFormats.MAPPING_OUTPUT.format(id="{identification}", mapper="bowtie")
    run:
        bowtie_mapper = MapperFactory.create(config[PipelineKeys.MAPPING][wildcards.identification][MappingKeys.LIBRARY])
        bowtie_mapper.map(config[PipelineKeys.MAPPING][wildcards.identification])
