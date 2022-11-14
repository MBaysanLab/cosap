from typing import List, Dict
import os
from cosap.variant_callers._variant_factory import VariantCallerFactory
from cosap._formats import FolderedOutputs
from cosap._pipeline_config import (
    VariantCallingKeys,
    PipelineKeys,
    SnakemakeConstraints,
)


def get_bams(wildcards):
    input_bams = []
    if (
        VariantCallingKeys.GERMLINE_INPUT
        in config[PipelineKeys.VARIANT_CALLING][wildcards.identification].keys()
    ):
        input_bams.append(
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
                VariantCallingKeys.GERMLINE_INPUT
            ]
        )

    if (
        VariantCallingKeys.TUMOR_INPUT
        in config[PipelineKeys.VARIANT_CALLING][wildcards.identification].keys()
    ):
        input_bams.append(
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
                VariantCallingKeys.TUMOR_INPUT
            ]
        )

    return input_bams


ruleorder: py2_variant_caller > variant_caller


rule variant_caller:
    input:
        bams=get_bams,
    output:
        vcf=FolderedOutputs.VARIANT_CALLING_OUTPUT,
    run:
        variant_caller = VariantCallerFactory.create(
            caller_type=config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
                VariantCallingKeys.LIBRARY
            ]
        )
        variant_caller.call_variants(
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification]
        )


rule variant_caller_with_gvcf_output:
    input:
        bams=get_bams,
    output:
        gvcf=FileFormats.GVCF_OUTPUT,
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
        vcf=FolderedOutputs.VARIANT_CALLING_OUTPUT,
    conda:
        "../../environments/py2_environment.yaml"
    wildcard_constraints:
        identification=SnakemakeConstraints.PY2_VARIANT_CALLERS,
    script:
        "../scripts/_py2_variantcaller.py"
