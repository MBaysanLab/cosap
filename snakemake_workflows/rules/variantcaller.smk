from typing import List, Dict
import os
from cosap.variant_callers._variant_factory import VariantCallerFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import VariantCallingKeys, PipelineKeys, SnakemakeConstraints

ruleorder: py2_variant_caller > variant_caller

rule variant_caller:
    input:
        germline_bam=lambda wildcards: config[PipelineKeys.VARIANT_CALLING][
            wildcards.identification
        ][VariantCallingKeys.GERMLINE_INPUT],
        tumor_bam=lambda wildcards: config[PipelineKeys.VARIANT_CALLING][
            wildcards.identification
        ][VariantCallingKeys.TUMOR_INPUT],
    output:
        vcf=FileFormats.GATK_SNP_OUTPUT,
    resources: variant_caller=1
    run:
        variant_caller = VariantCallerFactory.create(
            caller_type=config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
                VariantCallingKeys.LIBRARY
            ]
        )
        variant_caller.call_variants(
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification]
        )


rule py2_variant_caller:
    input:
        germline_bam=lambda wildcards: config[PipelineKeys.VARIANT_CALLING][
            wildcards.identification
        ][VariantCallingKeys.GERMLINE_INPUT],
        tumor_bam=lambda wildcards: config[PipelineKeys.VARIANT_CALLING][
            wildcards.identification
        ][VariantCallingKeys.TUMOR_INPUT],
    output:
        vcf=FileFormats.GATK_SNP_OUTPUT,
    resources: variant_caller=1
    conda: 
        "../../environments/py2_environment.yaml"
    wildcard_constraints:
        identification = SnakemakeConstraints.PY2_VARIANT_CALLERS
    script:
        "../scripts/_py2_variantcaller.py"