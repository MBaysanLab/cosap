from unittest.mock import Mock

from cosap._pipeline_config import BaseRecalibratorKeys, MDUPKeys, VariantCallingKeys
from cosap.tools.preprocessors._base_recalibrator import BaseRecalibrator
from cosap.tools.preprocessors._mark_duplicate import MarkDuplicate
from cosap.tools.variant_callers._haplotypecaller_variantcaller import (
    HaplotypeCallerVariantCaller,
)
from cosap.tools.variant_callers._mutect2_variantcaller import Mutect2VariantCaller


def test_markduplicate_parabricks_outputs_use_output_dir():
    config = {
        MDUPKeys.INPUT: "BAM/bwa/unprocessed.bam",
        MDUPKeys.OUTPUT: "mdup_sample.bam",
        MDUPKeys.OUTPUT_DIR: "PREPROCESSOR/mdup",
    }

    output_bam = MarkDuplicate._create_parabricks_output_path(config)
    queryname_bam = MarkDuplicate._create_parabricks_intermediate_path(
        config, "queryname"
    )
    markdup_bam = MarkDuplicate._create_parabricks_intermediate_path(config, "markdup")
    sort_command = MarkDuplicate.create_parabricks_sort_command(
        config[MDUPKeys.INPUT],
        queryname_bam,
        Mock(REF_FASTA="ref.fasta"),
        "queryname",
    )
    markdup_command = MarkDuplicate.create_parabricks_markdup_command(
        library_paths=Mock(REF_FASTA="ref.fasta"),
        app_config=Mock(),
        mdup_config=config,
        memory_handler=None,
        input_bam=queryname_bam,
        output_bam=markdup_bam,
        metrics_path=f"{output_bam}.duplicate_metrics",
    )

    assert output_bam == "PREPROCESSOR/mdup/mdup_sample.bam"
    assert sort_command[sort_command.index("--out-bam") + 1] == queryname_bam
    assert markdup_command[markdup_command.index("--in-bam") + 1] == queryname_bam
    assert markdup_command[markdup_command.index("--out-bam") + 1] == markdup_bam
    assert (
        markdup_command[markdup_command.index("--out-duplicate-metrics") + 1]
        == "PREPROCESSOR/mdup/mdup_sample.bam.duplicate_metrics"
    )


def test_base_recalibrator_parabricks_outputs_use_output_dir():
    config = {
        BaseRecalibratorKeys.INPUT: "PREPROCESSOR/mdup/mdup_sample.bam",
        BaseRecalibratorKeys.TABLE: "sample.table",
        BaseRecalibratorKeys.OUTPUT: "calibrated_sample.bam",
        BaseRecalibratorKeys.OUTPUT_DIR: "CALIBRATED_BAM",
    }

    bqsr_command = BaseRecalibrator._create_parabricks_bqsr_command(
        config,
        Mock(
            REF_FASTA="ref.fasta",
            MILLS_INDEL="mills.vcf",
            DBSNP="dbsnp.vcf",
            ONE_THOUSAND_G="1000g.vcf",
        ),
    )
    apply_command = BaseRecalibrator._create_parabricks_applybqsr_command(
        config, Mock(REF_FASTA="ref.fasta")
    )

    assert bqsr_command[bqsr_command.index("--out-recal-file") + 1] == (
        "CALIBRATED_BAM/sample.table"
    )
    assert apply_command[apply_command.index("--in-recal-file") + 1] == (
        "CALIBRATED_BAM/sample.table"
    )
    assert apply_command[apply_command.index("--out-bam") + 1] == (
        "CALIBRATED_BAM/calibrated_sample.bam"
    )


def test_mutect2_parabricks_output_uses_output_dir():
    config = {
        VariantCallingKeys.TUMOR_INPUT: "PREPROCESSOR/mdup/tumor.bam",
        VariantCallingKeys.GERMLINE_INPUT: "PREPROCESSOR/mdup/normal.bam",
        VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: "all_unfiltered_sample.vcf.gz",
        VariantCallingKeys.OUTPUT_DIR: "VCF/mutect2",
        VariantCallingKeys.PARAMS: {
            VariantCallingKeys.TUMOR_SAMPLE_NAME: "tumor",
            VariantCallingKeys.GERMLINE_SAMPLE_NAME: "normal",
        },
    }

    command = Mutect2VariantCaller._create_parabricks_mutectcaller_command(
        config, Mock(REF_FASTA="ref.fasta")
    )

    assert command[command.index("--out-vcf") + 1] == (
        "VCF/mutect2/all_unfiltered_sample.vcf.gz"
    )


def test_haplotypecaller_parabricks_output_uses_output_dir():
    config = {
        VariantCallingKeys.GERMLINE_INPUT: "PREPROCESSOR/mdup/normal.bam",
        VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: "all_unfiltered_sample.vcf.gz",
        VariantCallingKeys.GVCF_OUTPUT: "sample.g.vcf.gz",
        VariantCallingKeys.OUTPUT_DIR: "VCF/haplotypecaller",
        VariantCallingKeys.OUTPUT_TYPE: "VCF",
    }

    command = HaplotypeCallerVariantCaller._create_parabricks_haplotypecaller_command(
        config, Mock(REF_FASTA="ref.fasta")
    )

    assert command[command.index("--out-variants") + 1] == (
        "VCF/haplotypecaller/all_unfiltered_sample.vcf.gz"
    )
