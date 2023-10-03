import os
from cosap.tools.variant_callers._variant_factory import VariantCallerFactory
from cosap.tools.gene_fusion_callers._gene_fusion_caller_factory import GeneFusionCallerFactory
from cosap.tools.msi_callers._msi_caller_factory import MSICallerFactory
from cosap.tools.cnv_callers._cnv_caller_factory import CNVCallerFactory
from cosap._formats import FolderedOutputs
from cosap._pipeline_config import (
    VariantCallingKeys,
    PipelineKeys,
    SnakemakeConstraints,
    GeneFusionCallingKeys,
    MSICallingKeys,
    CNVCallingKeys,
)

def get_bams(wildcards, calling_rule: str) -> list[str]:
    input_bams = []
    if calling_rule == "variant_caller":

        normal_bam = config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
                VariantCallingKeys.GERMLINE_INPUT
            ] if VariantCallingKeys.GERMLINE_INPUT in config[PipelineKeys.VARIANT_CALLING][wildcards.identification].keys() else None
        tumor_bam = config[PipelineKeys.VARIANT_CALLING][wildcards.identification][
                VariantCallingKeys.TUMOR_INPUT
            ] if VariantCallingKeys.TUMOR_INPUT in config[PipelineKeys.VARIANT_CALLING][wildcards.identification].keys() else None

    elif calling_rule == "msi_caller":

        normal_bam = config[PipelineKeys.MSI][wildcards.identification][
                MSICallingKeys.NORMAL_INPUT
            ] if MSICallingKeys.NORMAL_INPUT in config[PipelineKeys.MSI][wildcards.identification].keys() else None
        tumor_bam = config[PipelineKeys.MSI][wildcards.identification][
                MSICallingKeys.TUMOR_INPUT
            ] if MSICallingKeys.TUMOR_INPUT in config[PipelineKeys.MSI][wildcards.identification].keys() else None
        
    elif calling_rule == "cnv_caller":

        normal_bam = config[PipelineKeys.CNV][wildcards.identification][
                CNVCallingKeys.NORMAL_INPUT
            ] if CNVCallingKeys.NORMAL_INPUT in config[PipelineKeys.CNV][wildcards.identification].keys() else None
        tumor_bam = config[PipelineKeys.CNV][wildcards.identification][
                CNVCallingKeys.TUMOR_INPUT
            ] if CNVCallingKeys.TUMOR_INPUT in config[PipelineKeys.CNV][wildcards.identification].keys() else None

    if normal_bam:
        input_bams.append(normal_bam)
    
    if tumor_bam:
        input_bams.append(tumor_bam)

    return input_bams


ruleorder: py2_variant_caller > variant_caller


rule variant_caller:
    input:
        lambda wildcards: get_bams(wildcards,"variant_caller"),
    output:
        vcf=FolderedOutputs.VARIANT_CALLING_VCF_OUTPUT,
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
        lambda wildcards: get_bams(wildcards,"variant_caller"),
    output:
        gvcf=FolderedOutputs.VARIANT_CALLING_GVCF_OUTPUT,
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
        lambda wildcards: get_bams(wildcards,"variant_caller"),
    output:
        vcf=FolderedOutputs.VARIANT_CALLING_VCF_OUTPUT,
    conda:
        "../../environments/py2_environment.yaml"
    wildcard_constraints:
        identification=SnakemakeConstraints.PY2_VARIANT_CALLERS,
    script:
        "../scripts/variantcaller_script.py"

rule genefusion_caller:
    input:
        fastqfiles=lambda wildcards: expand(
            config[PipelineKeys.GENEFUSION][wildcards.identification][
                GeneFusionCallingKeys.INPUT
            ].values()
        ),
    output:
        vcf=FolderedOutputs.GENEFUSION_OUTPUT,
    run:
        genefusion_caller = GeneFusionCallerFactory.create(
            caller_type=config[PipelineKeys.GENEFUSION][wildcards.identification][
                GeneFusionCallingKeys.LIBRARY
            ]
        )
        genefusion_caller.call(
            config[PipelineKeys.GENEFUSION][wildcards.identification]
        )

rule msi_caller:
    input:
        lambda wildcards: get_bams(wildcards,"msi_caller"),
    output:
        vcf=FolderedOutputs.MSI_OUTPUT,
    run:
        msi_caller = MSICallerFactory.create(
            caller_type=config[PipelineKeys.MSI][wildcards.identification][
                MSICallingKeys.LIBRARY
            ]
        )
        msi_caller.call(
            config[PipelineKeys.MSI][wildcards.identification]
        )

rule cnv_caller:
    input:
        lambda wildcards: get_bams(wildcards,"cnv_caller"),
    output:
        vcf=FolderedOutputs.CNV_OUTPUT,
    run:
        cnv_caller = CNVCallerFactory.create(
            caller_type=config[PipelineKeys.CNV][wildcards.identification][
                CNVCallingKeys.LIBRARY
            ]
        )
        cnv_caller.call(
            config[PipelineKeys.CNV][wildcards.identification]
        )