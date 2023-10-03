import os
from cosap.tools.preprocessors import PreprocessorFactory
from cosap.tools.quality_controllers import QualityControllerFactory
from cosap._formats import FileFormats, FolderedOutputs
from cosap._pipeline_config import (
    MDUPKeys,
    BaseRecalibratorKeys,
    TrimmingKeys,
    PipelineKeys,
    ElprepKeys,
    QualityControlKeys,
    IndexingKeys,
    SortingKeys,
)
from collections import defaultdict


rule fastp_trim:
    output:
        fastq_outputs=expand(
            FolderedOutputs.TRIMMING_OUTPUT,
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
        mdup_bam=FolderedOutputs.MDUP_OUTPUT,
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
        calibrated_bam=FolderedOutputs.CALIBRATION_OUTPUT,
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
        calibrated_bam=FolderedOutputs.ELPREP_CALIBRATION_OUTPUT,
    run:
        duplicate_remover = PreprocessorFactory.create(preprocessor_type="elprep")
        duplicate_remover.run_preprocessor(
            config[PipelineKeys.ELPREP_PROCESS][wildcards.identification]
        )


rule quality_control:
    input:
        bam=lambda wildcards: config[PipelineKeys.QUALITY_CONTROL][
            wildcards.identification
        ][QualityControlKeys.INPUT],
    output:
        qc=directory(FolderedOutputs.BAMQC_OUTPUT),
    run:
        quality_controller = QualityContollerFactory.create(
            quality_controller_type=config[PipelineKeys.QUALITY_CONTROL][
                wildcards.identification
            ][QualityControlKeys.LIBRARY]
        )
        quality_controller.run_qualitycontroller(
            config[PipelineKeys.QUALITY_CONTROL][wildcards.identification]
        )


rule bam_sorting:
    input:
        bam=lambda wildcards: config[PipelineKeys.SORTING][wildcards.identification][
            SortingKeys.INPUT
        ],
    output:
        sorted_bam=FileFormats.SORTING_OUTPUT,
    run:
        sorter = PreprocessorFactory.create(preprocessor_type="sorter")
        sorter.run_preprocessor(config[PipelineKeys.SORTING][wildcards.identification])


rule bam_indexing:
    input:
        bam=lambda wildcards: config[PipelineKeys.INDEX][wildcards.bam_file][
            IndexingKeys.INPUT
        ],
    output:
        index=FileFormats.INDEXING_OUTPUT,
    run:
        indexer = PreprocessorFactory.create(preprocessor_type="indexer")
        indexer.run_preprocessor(config[PipelineKeys.INDEX][wildcards.bam_file])
