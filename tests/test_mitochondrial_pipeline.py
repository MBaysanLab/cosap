import pytest

pytest.importorskip("pysam")

from cosap._pipeline_config import MitochondrialCallingKeys, PipelineKeys
from cosap.pipeline_builder.builders import MitochondrialVariantCaller
from cosap.workflows import MitochondrialPipelineInput


class _BamStep:
    name = "sample"
    next_step = None

    @staticmethod
    def get_output():
        return "CALIBRATED_BAM/calibrated_sample.bam"


def test_mitochondrial_builder_creates_dedicated_config():
    input_step = _BamStep()
    caller = MitochondrialVariantCaller(
        input_bam=input_step,
        input_bam_index="sample.bam.bai",
        median_autosomal_coverage=31.5,
        min_allele_fraction=0.01,
    )

    config = caller.get_config()
    core = config[PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING]["sample_mitochondrial"]

    assert input_step.next_step is caller
    assert core[MitochondrialCallingKeys.INPUT_BAM] == input_step.get_output()
    assert core[MitochondrialCallingKeys.INPUT_BAM_INDEX] == "sample.bam.bai"
    assert core[MitochondrialCallingKeys.CONTIG_NAME] == "chrM"
    assert core[MitochondrialCallingKeys.FILTERED_VCF].endswith(".chrM.vcf.gz")


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"INPUT_BAM": "sample.bam", "FASTQS": ("r1.fq.gz", "r2.fq.gz")},
        {"INPUT_BAM": "sample.cram", "INPUT_BAM_INDEX": "sample.crai"},
        {"INPUT_BAM": "sample.bam"},
        {"FASTQS": ("r1.fq.gz",)},
        {"FASTQS": ("r1.fq.gz", None)},
    ],
)
def test_mitochondrial_input_rejects_invalid_input_modes(kwargs):
    with pytest.raises(ValueError):
        MitochondrialPipelineInput(WORKDIR=".", SAMPLE_NAME="sample", **kwargs)


def test_mitochondrial_input_accepts_paired_fastqs():
    workflow_input = MitochondrialPipelineInput(
        WORKDIR=".", SAMPLE_NAME="sample", FASTQS=("r1.fq.gz", "r2.fq.gz")
    )

    assert workflow_input.MAPPER == "bwa"
    assert not workflow_input.SKIP_RECALIBRATION


def test_mitochondrial_builder_rejects_negative_f_score_beta():
    with pytest.raises(ValueError, match="f_score_beta"):
        MitochondrialVariantCaller(input_bam=_BamStep(), f_score_beta=-0.1)
