from typing import List, Dict
import os
from cosap.variant_callers._variant_factory import VariantCallerFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import VariantCallingKeys, PipelineKeys


rule variant_caller:
    input:
        germline_bam=lambda wildcards: config[PipelineKeys.VARIANT_CALLING][wildcards.identification][VariantCallingKeys.GERMLINE_INPUT],
        tumor_bam=lambda wildcards: config[PipelineKeys.VARIANT_CALLING][wildcards.identification][VariantCallingKeys.TUMOR_INPUT]
    output:
        vcf=config[PipelineKeys.VARIANT_CALLING][VariantCallingKeys.SNAKEMAKE_OUTPUT]
    run:
        variant_caller = VariantCallerFactory.create(
            caller_type=config[PipelineKeys.VARIANT_CALLING][wildcards.identification][VariantCallingKeys.LIBRARY]
        )
        variant_caller.call_variants(
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification]
        )
