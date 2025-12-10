"""
Example: Running GATK GermlineCNVCaller Pipeline

This example demonstrates how to use the GATK GermlineCNVCaller workflow
in COSAP for detecting rare germline copy number variants from multiple samples.

The GermlineCNVCaller pipeline implements the complete GATK workflow:
Phase 1 - Cohort Model Building:
    1. PreprocessIntervals - Prepare intervals (bin WGS or pad exomes)
    2. AnnotateIntervals - Annotate with GC content
    3. CollectReadCounts - Collect per-sample read counts  
    4. FilterIntervals - Filter problematic intervals
    5. DetermineGermlineContigPloidy - Determine contig ploidies

Phase 2 - Sample CNV Calling:
    6. GermlineCNVCaller - Call CNVs
    7. PostprocessGermlineCNVCalls - Generate per-sample VCFs

References:
- GATK Documentation: https://gatk.broadinstitute.org/hc/en-us/articles/360035531152
- Whitepaper: https://github.com/broadinstitute/gatk/blob/master/docs/CNV/germline-cnv-caller-model.pdf
"""

from cosap import (
    MDUP,
    FastqReader,
    Mapper,
    Pipeline,
    PipelineRunner,
    Recalibrator,
    Trimmer,
)
from cosap.pipeline_builder.builders import (
    GermlineCNVCaller,
    GermlineCNVCallerExome,
    GermlineCNVCallerSensitive,
)


def example_wgs_cohort():
    """
    Example: WGS cohort analysis with default parameters.
    
    This is the standard workflow for whole-genome sequencing data.
    Recommended minimum cohort size is 100 samples for optimal model training.
    """
    
    # Define samples (minimum 3 samples for cohort mode)
    samples = [
        {"name": "sample1", "r1": "./data/sample1_R1.fastq.gz", "r2": "./data/sample1_R2.fastq.gz"},
        {"name": "sample2", "r1": "./data/sample2_R1.fastq.gz", "r2": "./data/sample2_R2.fastq.gz"},
        {"name": "sample3", "r1": "./data/sample3_R1.fastq.gz", "r2": "./data/sample3_R2.fastq.gz"},
    ]
    
    # Create preprocessing steps for each sample
    calibrated_bams = []
    pipeline = Pipeline()
    
    for sample in samples:
        # Read FASTQ files
        fastqs = [
            FastqReader(sample["r1"], name=sample["name"], read=1),
            FastqReader(sample["r2"], name=sample["name"], read=2),
        ]
        
        # Trim reads
        trimmer = Trimmer(input_step=fastqs, name=sample["name"])
        pipeline.add(trimmer)
        
        # Map to reference
        mapper = Mapper(
            library="bwa2",
            input_step=trimmer,
            params={
                "read_groups": {
                    "ID": sample["name"],
                    "SM": sample["name"],
                    "PU": "unit1",
                    "PL": "ILLUMINA",
                    "LB": "lib1"
                }
            }
        )
        pipeline.add(mapper)
        
        # Mark duplicates
        mdup = MDUP(input_step=mapper)
        pipeline.add(mdup)
        
        # Base quality score recalibration
        bqsr = Recalibrator(input_step=mdup)
        pipeline.add(bqsr)
        
        calibrated_bams.append(bqsr)
    
    # Create GermlineCNVCaller pipeline with all calibrated BAMs
    # WGS mode: default bin_length=1000, padding=0
    germline_cnv = GermlineCNVCaller(
        input_steps=calibrated_bams,
        name="wgs_cohort_cnv",
        contig_ploidy_priors="./reference/contig_ploidy_priors.tsv",
        # For WGS, default parameters are fine:
        # bin_length=1000 (1kb bins)
        # padding=0 (no padding)
    )
    pipeline.add(germline_cnv)
    
    # Build and run pipeline
    config = pipeline.build(workdir="./workdir")
    
    runner = PipelineRunner()
    runner.run_pipeline(config)


def example_wgs_sensitive():
    """
    Example: WGS analysis with increased sensitivity parameters.
    
    Use this for detecting smaller CNV events. May require additional
    post-calling filtering due to increased false positives.
    """
    
    # Assume we have preprocessed BAMs
    from cosap.pipeline_builder.builders import BamReader
    
    calibrated_bams = [
        BamReader("./data/sample1.calibrated.bam", name="sample1"),
        BamReader("./data/sample2.calibrated.bam", name="sample2"),
        BamReader("./data/sample3.calibrated.bam", name="sample3"),
    ]
    
    # Use the sensitive preset for WGS
    germline_cnv = GermlineCNVCallerSensitive(
        input_steps=calibrated_bams,
        name="wgs_sensitive_cnv",
        contig_ploidy_priors="./reference/contig_ploidy_priors.tsv",
    )
    
    pipeline = Pipeline()
    pipeline.add(germline_cnv)
    
    config = pipeline.build(workdir="./workdir")
    
    runner = PipelineRunner()
    runner.run_pipeline(config)


def example_exome_cohort():
    """
    Example: Targeted exome analysis.
    
    For exome/panel data:
    - Use bin_length=0 (no binning)
    - Use padding=250 (pad target regions)
    - Provide target intervals file
    - Keep enable_bias_factors=True
    """
    
    from cosap.pipeline_builder.builders import BamReader
    
    calibrated_bams = [
        BamReader("./data/sample1.calibrated.bam", name="sample1"),
        BamReader("./data/sample2.calibrated.bam", name="sample2"),
        BamReader("./data/sample3.calibrated.bam", name="sample3"),
    ]
    
    # Use the exome preset
    germline_cnv = GermlineCNVCallerExome(
        input_steps=calibrated_bams,
        intervals="./reference/exome_targets.interval_list",  # Required for exomes
        contig_ploidy_priors="./reference/contig_ploidy_priors.tsv",
        name="exome_cohort_cnv",
    )
    
    pipeline = Pipeline()
    pipeline.add(germline_cnv)
    
    config = pipeline.build(workdir="./workdir")
    
    runner = PipelineRunner()
    runner.run_pipeline(config)


def example_custom_parameters():
    """
    Example: Custom parameter tuning.
    
    Advanced users can fine-tune all parameters for their specific data.
    """
    
    from cosap.pipeline_builder.builders import BamReader
    
    calibrated_bams = [
        BamReader("./data/sample1.calibrated.bam", name="sample1"),
        BamReader("./data/sample2.calibrated.bam", name="sample2"),
        BamReader("./data/sample3.calibrated.bam", name="sample3"),
    ]
    
    # Full custom configuration
    germline_cnv = GermlineCNVCaller(
        input_steps=calibrated_bams,
        intervals="./reference/targets.interval_list",
        contig_ploidy_priors="./reference/contig_ploidy_priors.tsv",
        name="custom_cnv",
        bin_length=500,  # Smaller bins for higher resolution
        padding=100,
        allosomal_contigs=["chrX", "chrY"],  # Sex chromosomes
        params={
            # GC content filtering
            "minimum_gc_content": 0.1,
            "maximum_gc_content": 0.9,
            
            # Mappability filtering (requires mappability track)
            # "mappability_track": "./reference/mappability.bed.gz",
            # "minimum_mappability": 0.9,
            # "maximum_mappability": 1.0,
            
            # Low count filtering
            "low_count_filter_count_threshold": 10,
            "low_count_filter_percentage_of_samples": 50.0,
            
            # Extreme count filtering
            "extreme_count_filter_minimum_percentile": 1.0,
            "extreme_count_filter_maximum_percentile": 99.0,
            "extreme_count_filter_percentage_of_samples": 90.0,
            
            # CNV calling sensitivity parameters
            "class_coherence_length": 5000.0,  # Default is 10000
            "cnv_coherence_length": 5000.0,    # Default is 10000
            "enable_bias_factors": True,        # Keep True for exomes
            
            # Advanced tuning
            "max_copy_number": 5,
            "p_active": 0.01,
            "p_alt": 1e-6,
            
            # Verbosity
            "verbosity": "DEBUG",  # For detailed logging
        },
    )
    
    pipeline = Pipeline()
    pipeline.add(germline_cnv)
    
    config = pipeline.build(workdir="./workdir")
    
    runner = PipelineRunner()
    runner.run_pipeline(config)


if __name__ == "__main__":
    # Run the WGS cohort example
    example_wgs_cohort()
