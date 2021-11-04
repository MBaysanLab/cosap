import os
from cosap.preprocessors._preprocessor_factory import PreprocessorFactory
from cosap._formats import FileFormats
from cosap._pipeline_config import (
    MDUPKeys,
    BaseRecalibratorKeys,
    TrimmingKeys,
    PipelineKeys,
)
from collections import defaultdict


rule fastp_trim:
    output:
        fastq_outputs=expand(
            config[PipelineKeys.TRIM][TrimmingKeys.SNAKEMAKE_OUTPUT], pair=["1", "2"]
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
        mdup_bam=config[PipelineKeys.MDUP][MDUPKeys.SNAKEMAKE_OUTPUT],
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
        calibrated_bam=config[PipelineKeys.CALIBRATE][
            BaseRecalibratorKeys.SNAKEMAKE_OUTPUT
        ],
    run:
        duplicate_remover = PreprocessorFactory.create(
            preprocessor_type="base_recalibrator"
        )
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.CALIBRATE][wildcards.identification]
        )
