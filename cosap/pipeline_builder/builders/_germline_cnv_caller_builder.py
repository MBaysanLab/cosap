"""
Pipeline Builder for GATK GermlineCNVCaller

This module provides the GermlineCNVCaller builder class for configuring
the GATK GermlineCNVCaller pipeline in COSAP.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import GermlineCNVCallingKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class GermlineCNVCaller(_IPipelineStep, _PipelineStep):
    """
    Pipeline builder for GATK GermlineCNVCaller workflow.
    
    This builder configures the complete GATK GermlineCNVCaller pipeline
    for detecting rare germline copy number variants from multiple samples.
    
    The pipeline implements the following workflow (cohort mode):
    Phase 1 - Cohort Model Building:
        1. PreprocessIntervals - Prepare intervals (bin WGS or pad exomes)
        2. AnnotateIntervals - Annotate with GC content (optional but recommended)
        3. CollectReadCounts - Collect per-sample read counts
        4. FilterIntervals - Filter problematic intervals
        5. DetermineGermlineContigPloidy - Determine contig ploidies
    
    Phase 2 - Sample CNV Calling:
        6. GermlineCNVCaller - Call CNVs
        7. PostprocessGermlineCNVCalls - Generate per-sample VCFs
    
    Args:
        input_steps: List of pipeline steps providing BAM files (e.g., calibrated BAMs)
        intervals: Path to interval file (BED or interval_list) for targeted analysis
        name: Optional name for this pipeline step
        contig_ploidy_priors: Path to contig ploidy priors file (required for cohort mode)
        run_mode: COHORT (default) or CASE
        bin_length: Bin size for WGS (default 1000bp, set to 0 for exomes)
        padding: Padding for exome intervals (default 0 for WGS, 250 for exomes)
        allosomal_contigs: List of allosomal contig names (default ["chrX", "chrY"])
        params: Additional parameters for fine-tuning
    
    Example:
        >>> bqsr1 = Recalibrator(input_step=mdup1)
        >>> bqsr2 = Recalibrator(input_step=mdup2)
        >>> bqsr3 = Recalibrator(input_step=mdup3)
        >>> 
        >>> germline_cnv = GermlineCNVCaller(
        ...     input_steps=[bqsr1, bqsr2, bqsr3],
        ...     intervals="targets.interval_list",
        ...     contig_ploidy_priors="contig_ploidy_priors.tsv",
        ... )
    """
    
    input_steps: List[_PipelineStep]
    intervals: str = None
    name: str = None
    contig_ploidy_priors: str = None
    run_mode: str = "COHORT"
    bin_length: int = 1000  # For WGS; set to 0 for exomes
    padding: int = 0  # For WGS; set to 250 for exomes
    allosomal_contigs: List[str] = field(default_factory=lambda: ["chrX", "chrY"])
    params: Dict = field(default_factory=dict)
    key: str = PipelineKeys.GERMLINE_CNV
    next_step: _PipelineStep = None
    output_dir: str = None
    library: str = "gatk_germline_cnv"

    def __post_init__(self):
        if self.name is None:
            # Generate name from input sample names
            name_parts = [step.name for step in self.input_steps[:3]]  # Limit to first 3
            if len(self.input_steps) > 3:
                name_parts.append(f"plus{len(self.input_steps) - 3}")
            self.name = "_".join(name_parts) + "_gcnv"
        
        # Link input steps to this step
        for step in self.input_steps:
            step.next_step = self
        
        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.GERMLINE_CNV, self.library)

    def _create_config(self) -> Dict:
        """Create configuration dictionary for the GATK GermlineCNVCaller pipeline."""
        
        # Collect input BAMs
        sample_inputs = [step.get_output() for step in self.input_steps]
        
        # Define output paths
        genotyped_intervals = FileFormats.GCNV_GENOTYPED_INTERVALS.format(
            identification=self.name
        )
        genotyped_segments = FileFormats.GCNV_GENOTYPED_SEGMENTS.format(
            identification=self.name
        )
        
        # Build params dict with all parameters
        all_params = {
            "bin_length": self.bin_length,
            "padding": self.padding,
            "allosomal_contigs": self.allosomal_contigs,
            **self.params,
        }
        
        config = {
            self.name: {
                GermlineCNVCallingKeys.LIBRARY: self.library,
                GermlineCNVCallingKeys.SAMPLE_INPUTS: sample_inputs,
                GermlineCNVCallingKeys.OUTPUT_DIR: self.output_dir,
                GermlineCNVCallingKeys.RUN_MODE: self.run_mode,
                GermlineCNVCallingKeys.GENOTYPED_INTERVALS: join_paths(
                    self.output_dir, "results", genotyped_intervals
                ),
                GermlineCNVCallingKeys.GENOTYPED_SEGMENTS: join_paths(
                    self.output_dir, "results", genotyped_segments
                ),
                GermlineCNVCallingKeys.PARAMS: all_params,
                GermlineCNVCallingKeys.LOG_FILE: self.log_file,
            }
        }
        
        # Add optional intervals
        if self.intervals:
            config[self.name][GermlineCNVCallingKeys.INTERVALS] = self.intervals
        
        # Add contig ploidy priors (required for cohort mode)
        if self.contig_ploidy_priors:
            config[self.name][GermlineCNVCallingKeys.CONTIG_PLOIDY_PRIORS] = self.contig_ploidy_priors
        
        return config

    def get_output(self) -> str:
        """Return the path to the genotyped intervals VCF output."""
        config = self.get_config()
        return config[self.key][self.name][GermlineCNVCallingKeys.GENOTYPED_INTERVALS]

    def get_segments_output(self) -> str:
        """Return the path to the genotyped segments VCF output."""
        config = self.get_config()
        return config[self.key][self.name][GermlineCNVCallingKeys.GENOTYPED_SEGMENTS]

    def get_config(self, core_only: bool = False) -> Dict:
        """
        Get the configuration dictionary.
        
        Args:
            core_only: If True, return only the step configuration without the pipeline key wrapper.
        
        Returns:
            Configuration dictionary for this pipeline step.
        """
        germline_cnv_config = self._create_config()
        if core_only:
            return germline_cnv_config[self.name]
        config = {self.key: germline_cnv_config}
        return config


@dataclass
class GermlineCNVCallerSensitive(GermlineCNVCaller):
    """
    GermlineCNVCaller with increased sensitivity parameters for WGS data.
    
    This preset applies the sensitivity-increasing parameters recommended
    in the GATK documentation for WGS data. These parameters dramatically
    increase the sensitivity of calling but may require additional
    post-calling filtering.
    
    From GATK documentation:
    - class_coherence_length: 1000.0 (default 10000)
    - cnv_coherence_length: 1000.0 (default 10000)
    - enable_bias_factors: False (default True)
    - interval_psi_scale: 1.0E-6 (default 0.001)
    - log_mean_bias_standard_deviation: 0.01 (default 0.1)
    - sample_psi_scale: 1.0E-6 (default 0.0001)
    """
    
    def __post_init__(self):
        # Set sensitive parameters
        self.params = {
            "class_coherence_length": 1000.0,
            "cnv_coherence_length": 1000.0,
            "enable_bias_factors": False,
            "interval_psi_scale": 1.0e-6,
            "log_mean_bias_standard_deviation": 0.01,
            "sample_psi_scale": 1.0e-6,
            **self.params,
        }
        super().__post_init__()


@dataclass
class GermlineCNVCallerExome(GermlineCNVCaller):
    """
    GermlineCNVCaller configured for targeted exome data.
    
    This preset applies the recommended parameters for exome/targeted
    panel data:
    - bin_length: 0 (disable binning, use provided intervals directly)
    - padding: 250 (pad target regions by 250bp on each side)
    - enable_bias_factors: True (recommended for exomes)
    """
    
    bin_length: int = 0
    padding: int = 250
    
    def __post_init__(self):
        # Ensure bias factors are enabled for exomes
        if "enable_bias_factors" not in self.params:
            self.params["enable_bias_factors"] = True
        super().__post_init__()
