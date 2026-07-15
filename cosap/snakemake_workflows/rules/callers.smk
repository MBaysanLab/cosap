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
    GermlineCNVCallingKeys,
    MitochondrialCallingKeys,
)
from cosap.tools.mitochondrial_callers import MitochondrialCallerFactory

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
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification],config["device"]
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
            config[PipelineKeys.VARIANT_CALLING][wildcards.identification], config["device"]
        )

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


def get_germline_cnv_bams(wildcards) -> list[str]:
    """Get all BAM inputs for germline CNV calling"""
    return config[PipelineKeys.GERMLINE_CNV][wildcards.identification][
        GermlineCNVCallingKeys.SAMPLE_INPUTS
    ]


rule germline_cnv_caller:
    """
    GATK GermlineCNVCaller pipeline for detecting rare germline CNVs.
    
    This rule handles the complete GATK GermlineCNVCaller workflow:
    1. PreprocessIntervals - Prepare analysis intervals
    2. AnnotateIntervals - Annotate with GC content
    3. CollectReadCounts - Collect per-sample read counts
    4. FilterIntervals - Filter problematic intervals
    5. DetermineGermlineContigPloidy - Determine contig ploidies
    6. GermlineCNVCaller - Call CNVs
    7. PostprocessGermlineCNVCalls - Generate final VCFs
    """
    input:
        bams=lambda wildcards: get_germline_cnv_bams(wildcards),
    output:
        intervals=FolderedOutputs.GCNV_GENOTYPED_INTERVALS,
        segments=FolderedOutputs.GCNV_GENOTYPED_SEGMENTS,
    run:
        cnv_caller = CNVCallerFactory.create(
            caller_type=config[PipelineKeys.GERMLINE_CNV][wildcards.identification][
                GermlineCNVCallingKeys.LIBRARY
            ]
        )
        cnv_caller.call(
            config[PipelineKeys.GERMLINE_CNV][wildcards.identification]
        )


def get_mitochondrial_inputs(wildcards) -> list[str]:
    caller_config = config[PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING][
        wildcards.identification
    ]
    inputs = [caller_config[MitochondrialCallingKeys.INPUT_BAM]]
    if MitochondrialCallingKeys.INPUT_BAM_INDEX in caller_config:
        inputs.append(caller_config[MitochondrialCallingKeys.INPUT_BAM_INDEX])
    return inputs


rule mitochondrial_variant_caller:
    """GATK dual-reference mitochondrial SNV/indel discovery."""
    input:
        get_mitochondrial_inputs,
    output:
        vcf=FolderedOutputs.MITOCHONDRIAL_FILTERED_VCF,
    run:
        caller = MitochondrialCallerFactory.create(
            config[PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING][
                wildcards.identification
            ][MitochondrialCallingKeys.LIBRARY]
        )
        caller.call(
            config[PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING][
                wildcards.identification
            ]
        )
