"""Public COSAP workflow for GATK mitochondrial short-variant discovery."""

from dataclasses import dataclass
from typing import Optional, Tuple

from ..pipeline_builder import (
    BamReader,
    FastqReader,
    MDUP,
    Mapper,
    MitochondrialVariantCaller,
    Pipeline,
    Recalibrator,
    Trimmer,
)
from ..runners.pipeline_runner import PipelineRunner


@dataclass
class MitochondrialPipelineInput:
    WORKDIR: str
    SAMPLE_NAME: str
    INPUT_BAM: Optional[str] = None
    INPUT_BAM_INDEX: Optional[str] = None
    FASTQS: Optional[Tuple[str, str]] = None
    MAPPER: str = "bwa"
    SKIP_RECALIBRATION: bool = False
    CONTIG_NAME: str = "chrM"
    MEDIAN_AUTOSOMAL_COVERAGE: Optional[float] = None
    MIN_ALLELE_FRACTION: Optional[float] = None
    F_SCORE_BETA: Optional[float] = None
    CONTAMINATION_ESTIMATE: Optional[float] = None
    RUN_HAPLOCHECKER: bool = True
    DEVICE: str = "cpu"

    def __post_init__(self):
        has_bam = self.INPUT_BAM is not None
        has_fastqs = self.FASTQS is not None
        if has_bam == has_fastqs:
            raise ValueError("Provide exactly one of INPUT_BAM or FASTQS.")
        if has_bam:
            if not self.INPUT_BAM.endswith(".bam"):
                raise ValueError("INPUT_BAM must be a BAM file; CRAM is not yet supported.")
            if not self.INPUT_BAM_INDEX:
                raise ValueError("INPUT_BAM_INDEX is required with INPUT_BAM.")
        if has_fastqs and (len(self.FASTQS) != 2 or not all(self.FASTQS)):
            raise ValueError("FASTQS must contain exactly non-empty R1 and R2 paths.")
        if self.DEVICE.lower() != "cpu":
            raise ValueError("The mitochondrial workflow currently supports CPU execution only.")


class MitochondrialPipeline:
    """Build and run mitochondrial calling from paired FASTQs or an aligned BAM."""

    def __init__(self, mitochondrial_input: MitochondrialPipelineInput):
        self.input = mitochondrial_input
        self.pipeline = Pipeline()
        self._create_pipeline()
        self.config = self.pipeline.build(workdir=self.input.WORKDIR)

    def _caller(self, input_bam):
        caller = MitochondrialVariantCaller(
            input_bam=input_bam,
            input_bam_index=self.input.INPUT_BAM_INDEX,
            name=self.input.SAMPLE_NAME,
            contig_name=self.input.CONTIG_NAME,
            median_autosomal_coverage=self.input.MEDIAN_AUTOSOMAL_COVERAGE,
            min_allele_fraction=self.input.MIN_ALLELE_FRACTION,
            f_score_beta=self.input.F_SCORE_BETA,
            contamination_estimate=self.input.CONTAMINATION_ESTIMATE,
            run_haplochecker=self.input.RUN_HAPLOCHECKER,
        )
        self.pipeline.add(caller)

    def _create_pipeline(self):
        if self.input.INPUT_BAM:
            self._caller(BamReader(self.input.INPUT_BAM, name=self.input.SAMPLE_NAME))
            return

        readers = [
            FastqReader(self.input.FASTQS[0], read=1, name=self.input.SAMPLE_NAME),
            FastqReader(self.input.FASTQS[1], read=2, name=self.input.SAMPLE_NAME),
        ]
        trimmer = Trimmer(input_step=readers, name=self.input.SAMPLE_NAME)
        mapper = Mapper(
            library=self.input.MAPPER,
            input_step=trimmer,
            name=self.input.SAMPLE_NAME,
            params={
                "read_groups": {
                    "ID": "0",
                    "SM": self.input.SAMPLE_NAME,
                    "PU": "0",
                    "PL": "ILLUMINA",
                    "LB": "0",
                }
            },
        )
        mdup = MDUP(input_step=mapper, name=self.input.SAMPLE_NAME)
        self.pipeline.add(trimmer)
        self.pipeline.add(mapper)
        self.pipeline.add(mdup)
        final_bam = mdup
        if not self.input.SKIP_RECALIBRATION:
            recalibrator = Recalibrator(input_step=mdup, name=self.input.SAMPLE_NAME)
            self.pipeline.add(recalibrator)
            final_bam = recalibrator
        self._caller(final_bam)

    def run_pipeline(self):
        return PipelineRunner(device="cpu").run_pipeline(self.config)
