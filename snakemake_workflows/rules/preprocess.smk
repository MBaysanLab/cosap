import os
from cosap.preprocessors._preprocessor_factory import PreprocessorFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import MappingKeys, PipelineKeys


rule fastp_trim:
    output:
        fastq=FileFormats.TRIMMING_OUTPUT.format(identification="{identification}"),
    run:
        trimmer = PreprocessorFactory.create(
            preprocessor_type = "trimmer"
        )
        trimmer.run_preprocessor(config[PipelineKeys.TRIM][wildcards.identification])

rule mark_dup:
    input:
        bam=lambda wildcards: FileFormats.MAPPING_OUTPUT.format(
            identification=wildcards.identification
        ),
    output:
        mdup_bam=FileFormats.MDUP_OUTPUT.format(identification="{identification}"),
    run:
        duplicate_remover = PreprocessorFactory.create(
            preprocessor_type="mark_duplicate"
        )
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.MDUP][wildcards.identification]
        )


rule gatk_base_cal:
    input:
        bam=lambda wildcards: FileFormats.MDUP_OUTPUT.format(
            identification=wildcards.identification
        ),
    output:
        calibrated_bam=FileFormats.CALIBRATION_OUTPUT.format(
            identification="{identification}"
        ),
    run:
        duplicate_remover = PreprocessorFactory.create(
            preprocessor_type="base_recalibrator"
        )
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.CALIBRATE][wildcards.identification]
        )
