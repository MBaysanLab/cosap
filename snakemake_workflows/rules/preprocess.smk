import os
from cosap.preprocessors._preprocessor_factory import PreprocessorFactory
from cosap._formats import SnakemakeOutputFormats
from cosap._pipeline_config import (
    MDUPKeys,
    BaseRecalibratorKeys,
    TrimmingKeys,
    PipelineKeys,
    ElprepKeys,
)
from collections import defaultdict


rule fastp_trim:
    output:
        fastq_outputs=expand(
            SnakemakeOutputFormats.TRIMMING_OUTPUT,
            pair=["1", "2"],
            identification="{identification}",
        ),
    run:
        trimmer = PreprocessorFactory.create(preprocessor_type="trimmer")
        trimmer.run_preprocessor(config[PipelineKeys.TRIM][wildcards.identification])


rule mark_dup:
    input:
        bam=lambda wildcards: config[PipelineKeys.MDUP][wildcards.identification][
            MDUPKeys.INPUT
        ],
    output:
        mdup_bam=SnakemakeOutputFormats.MDUP_OUTPUT,
    run:
        duplicate_remover = PreprocessorFactory.create(
            preprocessor_type="mark_duplicate"
        )
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.MDUP][wildcards.identification]
        )


rule gatk_base_cal:
    input:
        bam=lambda wildcards: config[PipelineKeys.CALIBRATE][wildcards.identification][
            BaseRecalibratorKeys.INPUT
        ],
    output:
        calibrated_bam=SnakemakeOutputFormats.CALIBRATION_OUTPUT,
    run:
        duplicate_remover = PreprocessorFactory.create(
            preprocessor_type="base_recalibrator"
        )
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.CALIBRATE][wildcards.identification]
        )


rule elprep_cal:
    input:
        bam=lambda wildcards: config[PipelineKeys.ELPREP_PROCESS][
            wildcards.identification
        ][ElprepKeys.INPUT],
    output:
        calibrated_bam=SnakemakeOutputFormats.ELPREP_CALIBRATION_OUTPUT,
    run:
        duplicate_remover = PreprocessorFactory.create(preprocessor_type="elprep")
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.ELPREP_PROCESS][wildcards.identification]
        )
