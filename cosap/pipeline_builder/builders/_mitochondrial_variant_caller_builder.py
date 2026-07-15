"""Builder for the GATK mitochondrial short-variant discovery workflow."""

from dataclasses import dataclass, field
from typing import Dict, Optional

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import MitochondrialCallingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class MitochondrialVariantCaller(_IPipelineStep, _PipelineStep):
    """Configure GATK's dual-reference mitochondrial Mutect2 workflow."""

    input_bam: _PipelineStep
    input_bam_index: Optional[str] = None
    name: Optional[str] = None
    contig_name: str = "chrM"
    median_autosomal_coverage: Optional[float] = None
    min_allele_fraction: Optional[float] = None
    f_score_beta: Optional[float] = None
    contamination_estimate: Optional[float] = None
    run_haplochecker: bool = True
    params: Dict = field(default_factory=dict)
    key: str = PipelineKeys.MITOCHONDRIAL_VARIANT_CALLING
    next_step: _PipelineStep = None
    output_dir: Optional[str] = None
    library: str = "gatk_mitochondrial"

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.input_bam.name}_mitochondrial"
        if self.input_bam_index is not None and not self.input_bam_index.endswith(
            (".bai", ".csi")
        ):
            raise ValueError("input_bam_index must be a BAM .bai or .csi index file.")
        if not self.contig_name:
            raise ValueError("contig_name is required.")
        if self.median_autosomal_coverage is not None and self.median_autosomal_coverage < 0:
            raise ValueError("median_autosomal_coverage must be non-negative.")
        if self.min_allele_fraction is not None and not 0 <= self.min_allele_fraction <= 1:
            raise ValueError("min_allele_fraction must be between 0 and 1.")
        if self.f_score_beta is not None and self.f_score_beta < 0:
            raise ValueError("f_score_beta must be non-negative.")
        if self.contamination_estimate is not None and not 0 <= self.contamination_estimate <= 1:
            raise ValueError("contamination_estimate must be between 0 and 1.")

        self.input_bam.next_step = self
        if self.output_dir is None:
            self.output_dir = join_paths(
                OutputFolders.MITOCHONDRIAL_VARIANT_CALLING, self.library
            )

    def _output(self, template: str) -> str:
        return template.format(identification=self.name)

    def _create_config(self) -> Dict:
        config = {
            MitochondrialCallingKeys.LIBRARY: self.library,
            MitochondrialCallingKeys.INPUT_BAM: self.input_bam.get_output(),
            MitochondrialCallingKeys.OUTPUT_DIR: self.output_dir,
            MitochondrialCallingKeys.CONTIG_NAME: self.contig_name,
            MitochondrialCallingKeys.RUN_HAPLOCHECKER: self.run_haplochecker,
            MitochondrialCallingKeys.PARAMS: self.params,
            MitochondrialCallingKeys.LOG_FILE: self.log_file,
            MitochondrialCallingKeys.SUBSET_BAM: self._output(FileFormats.MT_SUBSET_BAM),
            MitochondrialCallingKeys.UNMAPPED_BAM: self._output(FileFormats.MT_UNMAPPED_BAM),
            MitochondrialCallingKeys.STANDARD_ALIGNED_BAM: self._output(
                FileFormats.MT_STANDARD_ALIGNED_BAM
            ),
            MitochondrialCallingKeys.SHIFTED_ALIGNED_BAM: self._output(
                FileFormats.MT_SHIFTED_ALIGNED_BAM
            ),
            MitochondrialCallingKeys.STANDARD_RAW_VCF: self._output(
                FileFormats.MT_STANDARD_RAW_VCF
            ),
            MitochondrialCallingKeys.SHIFTED_RAW_VCF: self._output(
                FileFormats.MT_SHIFTED_RAW_VCF
            ),
            MitochondrialCallingKeys.MERGED_RAW_VCF: self._output(
                FileFormats.MT_MERGED_RAW_VCF
            ),
            MitochondrialCallingKeys.MERGED_STATS: self._output(
                FileFormats.MT_MERGED_STATS
            ),
            MitochondrialCallingKeys.FILTERED_VCF: self._output(
                FileFormats.MT_FILTERED_VCF
            ),
            MitochondrialCallingKeys.SPLIT_VCF: self._output(FileFormats.MT_SPLIT_VCF),
            MitochondrialCallingKeys.COVERAGE_METRICS: self._output(
                FileFormats.MT_COVERAGE_METRICS
            ),
            MitochondrialCallingKeys.DUPLICATE_METRICS: self._output(
                FileFormats.MT_DUPLICATE_METRICS
            ),
            MitochondrialCallingKeys.CONTAMINATION_REPORT: self._output(
                FileFormats.MT_CONTAMINATION_REPORT
            ),
        }
        if self.input_bam_index is not None:
            config[MitochondrialCallingKeys.INPUT_BAM_INDEX] = self.input_bam_index
        if self.median_autosomal_coverage is not None:
            config[MitochondrialCallingKeys.MEDIAN_AUTOSOMAL_COVERAGE] = (
                self.median_autosomal_coverage
            )
        if self.min_allele_fraction is not None:
            config[MitochondrialCallingKeys.MIN_ALLELE_FRACTION] = self.min_allele_fraction
        if self.f_score_beta is not None:
            config[MitochondrialCallingKeys.F_SCORE_BETA] = self.f_score_beta
        if self.contamination_estimate is not None:
            config[MitochondrialCallingKeys.CONTAMINATION_ESTIMATE] = (
                self.contamination_estimate
            )
        return config

    def get_config(self, core_only: bool = False) -> Dict:
        config = self._create_config()
        if core_only:
            return config
        return {self.key: {self.name: config}}

    def get_output(self) -> str:
        return join_paths(self.output_dir, self._output(FileFormats.MT_FILTERED_VCF))
