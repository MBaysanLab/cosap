from typing import List, Dict
import os
from cosap.annotators import AnnotatorFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import AnnotatorKeys, PipelineKeys


rule annotator:
    input:
        vcf=lambda wildcards: 
            config[PipelineKeys.ANNOTATION][wildcards.identification][AnnotatorKeys.INPUT]
    output:
        annotation=FileFormats.ANNOTATING_OUTPUT,
    run:
        annotator = AnnotatorFactory.create(
            config[PipelineKeys.ANNOTATION][wildcards.identification][AnnotatorKeys.LIBRARY]
        )
        annotator.annotate(config[PipelineKeys.ANNOTATION][wildcards.identification])
