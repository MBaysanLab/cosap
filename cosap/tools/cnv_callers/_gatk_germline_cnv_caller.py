"""
GATK GermlineCNVCaller Pipeline Implementation

This module implements the complete GATK GermlineCNVCaller workflow for detecting
rare germline copy number variants (CNVs) as described in:
https://gatk.broadinstitute.org/hc/en-us/articles/360035531152

The workflow consists of the following steps:
1. PreprocessIntervals - Prepare intervals for analysis (bin WGS or pad exomes)
2. AnnotateIntervals - Annotate intervals with GC content and mappability (optional but recommended)
3. CollectReadCounts - Collect read counts per interval per sample
4. FilterIntervals - Filter intervals based on annotations and extreme counts
5. DetermineGermlineContigPloidy - Determine baseline contig ploidies
6. GermlineCNVCaller - Call germline CNVs (cohort or case mode)
7. PostprocessGermlineCNVCalls - Generate per-sample VCF results
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import GermlineCNVCallingKeys
from ...runners.docker_runner import DockerRunner
from ._cnv_callers import _CNVCaller


class GATKGermlineCNVCaller(_CNVCaller):
    """
    GATK GermlineCNVCaller - Complete pipeline for germline CNV calling.
    
    This class orchestrates the entire GATK GermlineCNVCaller workflow,
    executing all steps in sequence to produce final VCF outputs.
    All GATK commands are executed inside Docker containers.
    """

    @classmethod
    def _create_preprocess_intervals_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """
        Create command for GATK PreprocessIntervals.
        
        For WGS: bins intervals with default 1000bp bin size
        For exomes: pads target regions with configurable padding (default 250bp)
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        bin_length = params.get("bin_length", 1000)  # Default for WGS
        padding = params.get("padding", 0)  # Default 0 for WGS, 250 for exomes
        
        command = [
            "gatk",
            "PreprocessIntervals",
            "-R", library_paths.REF_FASTA,
            "--bin-length", str(bin_length),
            "--padding", str(padding),
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "-O", caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS],
        ]
        
        # Add intervals if provided (required for exome, optional for WGS)
        if GermlineCNVCallingKeys.INTERVALS in caller_config:
            command.extend(["-L", caller_config[GermlineCNVCallingKeys.INTERVALS]])
        
        # Add exclude intervals if provided
        if params.get("exclude_intervals"):
            command.extend(["-XL", params["exclude_intervals"]])
        
        return command

    @classmethod
    def _create_annotate_intervals_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """
        Create command for GATK AnnotateIntervals.
        
        Annotates intervals with GC content and optionally mappability
        and segmental duplication content.
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        
        command = [
            "gatk",
            "AnnotateIntervals",
            "-R", library_paths.REF_FASTA,
            "-L", caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS],
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "-O", caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS],
        ]
        
        # Add optional mappability track
        if params.get("mappability_track"):
            command.extend(["--mappability-track", params["mappability_track"]])
        
        # Add optional segmental duplication track
        if params.get("segmental_duplication_track"):
            command.extend(["--segmental-duplication-track", params["segmental_duplication_track"]])
        
        return command

    @classmethod
    def _create_collect_read_counts_command(
        cls,
        bam_file: str,
        intervals: str,
        output_file: str,
        library_paths: LibraryPaths,
        params: Optional[Dict] = None,
    ) -> List[str]:
        """
        Create command for GATK CollectReadCounts.
        
        Collects read counts per interval for a single sample.
        Default format is HDF5 (more efficient for downstream tools).
        """
        params = params or {}
        
        command = [
            "gatk",
            "CollectReadCounts",
            "-I", bam_file,
            "-L", intervals,
            "-R", library_paths.REF_FASTA,
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "-O", output_file,
        ]
        
        # Optionally output TSV format instead of HDF5
        if params.get("format") == "TSV":
            command.extend(["--format", "TSV"])
        
        return command

    @classmethod
    def _create_filter_intervals_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """
        Create command for GATK FilterIntervals.
        
        Filters intervals based on GC content, mappability, and extreme counts.
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        
        command = [
            "gatk",
            "FilterIntervals",
            "-L", caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS],
            "--annotated-intervals", caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS],
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "-O", caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS],
        ]
        
        # Add all read count files
        read_counts = caller_config[GermlineCNVCallingKeys.READ_COUNTS]
        if isinstance(read_counts, dict):
            for rc_file in read_counts.values():
                command.extend(["-I", rc_file])
        else:
            for rc_file in read_counts:
                command.extend(["-I", rc_file])
        
        # GC content filtering thresholds
        if params.get("minimum_gc_content"):
            command.extend(["--minimum-gc-content", str(params["minimum_gc_content"])])
        if params.get("maximum_gc_content"):
            command.extend(["--maximum-gc-content", str(params["maximum_gc_content"])])
        
        # Mappability filtering thresholds
        if params.get("minimum_mappability"):
            command.extend(["--minimum-mappability", str(params["minimum_mappability"])])
        if params.get("maximum_mappability"):
            command.extend(["--maximum-mappability", str(params["maximum_mappability"])])
        
        # Low count filter thresholds
        if params.get("low_count_filter_count_threshold"):
            command.extend(["--low-count-filter-count-threshold", str(params["low_count_filter_count_threshold"])])
        if params.get("low_count_filter_percentage_of_samples"):
            command.extend(["--low-count-filter-percentage-of-samples", str(params["low_count_filter_percentage_of_samples"])])
        
        # Extreme count filter thresholds
        if params.get("extreme_count_filter_minimum_percentile"):
            command.extend(["--extreme-count-filter-minimum-percentile", str(params["extreme_count_filter_minimum_percentile"])])
        if params.get("extreme_count_filter_maximum_percentile"):
            command.extend(["--extreme-count-filter-maximum-percentile", str(params["extreme_count_filter_maximum_percentile"])])
        if params.get("extreme_count_filter_percentage_of_samples"):
            command.extend(["--extreme-count-filter-percentage-of-samples", str(params["extreme_count_filter_percentage_of_samples"])])
        
        # Exclude intervals
        if params.get("exclude_intervals"):
            command.extend(["-XL", params["exclude_intervals"]])
        
        return command

    @classmethod
    def _create_determine_ploidy_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """
        Create command for GATK DetermineGermlineContigPloidy.
        
        Determines baseline contig ploidies using sample coverages
        and contig ploidy priors.
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        run_mode = caller_config.get(GermlineCNVCallingKeys.RUN_MODE, "COHORT")
        
        command = [
            "gatk",
            "DetermineGermlineContigPloidy",
            "-L", caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS],
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--output", caller_config[GermlineCNVCallingKeys.OUTPUT_DIR],
            "--output-prefix", "ploidy",
        ]
        
        # Add read count files
        read_counts = caller_config[GermlineCNVCallingKeys.READ_COUNTS]
        if isinstance(read_counts, dict):
            for rc_file in read_counts.values():
                command.extend(["-I", rc_file])
        else:
            for rc_file in read_counts:
                command.extend(["-I", rc_file])
        
        # COHORT mode requires contig ploidy priors
        if run_mode == "COHORT":
            if GermlineCNVCallingKeys.CONTIG_PLOIDY_PRIORS in caller_config:
                command.extend(["--contig-ploidy-priors", caller_config[GermlineCNVCallingKeys.CONTIG_PLOIDY_PRIORS]])
            elif hasattr(library_paths, 'CONTIG_PLOIDY_PRIORS'):
                command.extend(["--contig-ploidy-priors", library_paths.CONTIG_PLOIDY_PRIORS])
        else:
            # CASE mode requires a model
            if params.get("ploidy_model"):
                command.extend(["--model", params["ploidy_model"]])
        
        # Optional exclude intervals
        if params.get("exclude_intervals"):
            command.extend(["-XL", params["exclude_intervals"]])
        
        # Verbosity for tracking model fitting
        command.extend(["--verbosity", params.get("verbosity", "INFO")])
        
        return command

    @classmethod
    def _create_germline_cnv_caller_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        scatter_index: Optional[int] = None,
    ) -> List[str]:
        """
        Create command for GATK GermlineCNVCaller.
        
        Calls germline CNVs using a denoising model. Can run in COHORT mode
        (builds model and calls simultaneously) or CASE mode (uses existing model).
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        run_mode = caller_config.get(GermlineCNVCallingKeys.RUN_MODE, "COHORT")
        
        output_prefix = "gcnv"
        if scatter_index is not None:
            output_prefix = f"gcnv_shard_{scatter_index}"
        
        command = [
            "gatk",
            "GermlineCNVCaller",
            "--run-mode", run_mode,
            "-L", caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS],
            "--contig-ploidy-calls", caller_config[GermlineCNVCallingKeys.PLOIDY_CALLS],
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--output", caller_config[GermlineCNVCallingKeys.OUTPUT_DIR],
            "--output-prefix", output_prefix,
        ]
        
        # Add read count files
        read_counts = caller_config[GermlineCNVCallingKeys.READ_COUNTS]
        if isinstance(read_counts, dict):
            for rc_file in read_counts.values():
                command.extend(["-I", rc_file])
        else:
            for rc_file in read_counts:
                command.extend(["-I", rc_file])
        
        # Add annotated intervals for explicit GC correction (recommended)
        if GermlineCNVCallingKeys.ANNOTATED_INTERVALS in caller_config:
            command.extend(["--annotated-intervals", caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS]])
        
        # CASE mode requires a model
        if run_mode == "CASE" and params.get("cnv_model"):
            command.extend(["--model", params["cnv_model"]])
        
        # Sensitivity tuning parameters for WGS (from GATK documentation)
        if params.get("class_coherence_length"):
            command.extend(["--class-coherence-length", str(params["class_coherence_length"])])
        if params.get("cnv_coherence_length"):
            command.extend(["--cnv-coherence-length", str(params["cnv_coherence_length"])])
        if params.get("enable_bias_factors") is not None:
            command.extend(["--enable-bias-factors", str(params["enable_bias_factors"]).lower()])
        if params.get("interval_psi_scale"):
            command.extend(["--interval-psi-scale", str(params["interval_psi_scale"])])
        if params.get("log_mean_bias_standard_deviation"):
            command.extend(["--log-mean-bias-standard-deviation", str(params["log_mean_bias_standard_deviation"])])
        if params.get("sample_psi_scale"):
            command.extend(["--sample-psi-scale", str(params["sample_psi_scale"])])
        if params.get("depth_correction_tau"):
            command.extend(["--depth-correction-tau", str(params["depth_correction_tau"])])
        if params.get("p_active"):
            command.extend(["--p-active", str(params["p_active"])])
        if params.get("p_alt"):
            command.extend(["--p-alt", str(params["p_alt"])])
        
        # Max copy number (default is 5)
        if params.get("max_copy_number"):
            command.extend(["--max-copy-number", str(params["max_copy_number"])])
        
        # Verbosity
        command.extend(["--verbosity", params.get("verbosity", "INFO")])
        
        return command

    @classmethod
    def _create_postprocess_calls_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        sample_index: int = 0,
        sample_name: Optional[str] = None,
    ) -> List[str]:
        """
        Create command for GATK PostprocessGermlineCNVCalls.
        
        Consolidates scattered GermlineCNVCaller results, performs segmentation,
        and generates per-sample VCF files.
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        
        # Handle output naming
        if sample_name:
            intervals_output = caller_config[GermlineCNVCallingKeys.GENOTYPED_INTERVALS].replace(
                "{sample}", sample_name
            )
            segments_output = caller_config[GermlineCNVCallingKeys.GENOTYPED_SEGMENTS].replace(
                "{sample}", sample_name
            )
            denoised_output = caller_config.get(
                GermlineCNVCallingKeys.DENOISED_COPY_RATIOS,
                intervals_output.replace("genotyped_intervals", "denoised_copy_ratios").replace(".vcf.gz", ".tsv")
            ).replace("{sample}", sample_name)
        else:
            intervals_output = caller_config[GermlineCNVCallingKeys.GENOTYPED_INTERVALS]
            segments_output = caller_config[GermlineCNVCallingKeys.GENOTYPED_SEGMENTS]
            denoised_output = caller_config.get(
                GermlineCNVCallingKeys.DENOISED_COPY_RATIOS,
                intervals_output.replace("genotyped_intervals", "denoised_copy_ratios").replace(".vcf.gz", ".tsv")
            )
        
        # Use REF_DICT if available, otherwise derive from REF_FASTA
        if hasattr(library_paths, 'REF_DICT'):
            sequence_dict = library_paths.REF_DICT
        else:
            sequence_dict = library_paths.REF_FASTA.replace(".fasta", ".dict").replace(".fa", ".dict")
        
        command = [
            "gatk",
            "PostprocessGermlineCNVCalls",
            "--sample-index", str(sample_index),
            "--output-genotyped-intervals", intervals_output,
            "--output-genotyped-segments", segments_output,
            "--output-denoised-copy-ratios", denoised_output,
            "--sequence-dictionary", sequence_dict,
        ]
        
        # Add model shard paths
        model_shards = caller_config.get(GermlineCNVCallingKeys.CNV_MODEL)
        if isinstance(model_shards, list):
            for shard in model_shards:
                command.extend(["--model-shard-path", shard])
        else:
            command.extend(["--model-shard-path", model_shards])
        
        # Add calls shard paths
        calls_shards = caller_config.get(GermlineCNVCallingKeys.CNV_CALLS)
        if isinstance(calls_shards, list):
            for shard in calls_shards:
                command.extend(["--calls-shard-path", shard])
        else:
            command.extend(["--calls-shard-path", calls_shards])
        
        # Add contig ploidy calls
        command.extend(["--contig-ploidy-calls", caller_config[GermlineCNVCallingKeys.PLOIDY_CALLS]])
        
        # Allosomal contigs (for proper handling of sex chromosomes)
        allosomal_contigs = params.get("allosomal_contigs", ["chrX", "chrY"])
        for contig in allosomal_contigs:
            command.extend(["--allosomal-contig", contig])
        
        # Autosomal reference copy number (default is 2)
        if params.get("autosomal_ref_copy_number"):
            command.extend(["--autosomal-ref-copy-number", str(params["autosomal_ref_copy_number"])])
        
        return command

    @classmethod
    def call(cls, caller_config: Dict, device: str = "cpu"):
        """
        Execute the complete GATK GermlineCNVCaller pipeline.
        
        Pipeline steps:
        1. PreprocessIntervals - Prepare analysis intervals
        2. AnnotateIntervals - Add GC content annotations
        3. CollectReadCounts - Collect per-sample read counts
        4. FilterIntervals - Filter problematic intervals
        5. DetermineGermlineContigPloidy - Determine contig ploidies
        6. GermlineCNVCaller - Call CNVs
        7. PostprocessGermlineCNVCalls - Generate final VCFs
        
        All commands are executed inside GATK Docker containers.
        Steps are skipped if their output files already exist.
        """
        library_paths = LibraryPaths()
        app_config = AppConfig()
        
        output_dir = caller_config[GermlineCNVCallingKeys.OUTPUT_DIR]
        os.makedirs(output_dir, exist_ok=True)
        
        # Create subdirectories
        intervals_dir = os.path.join(output_dir, "intervals")
        read_counts_dir = os.path.join(output_dir, "read_counts")
        ploidy_dir = os.path.join(output_dir, "ploidy")
        model_dir = os.path.join(output_dir, "model")
        calls_dir = os.path.join(output_dir, "calls")
        results_dir = os.path.join(output_dir, "results")
        
        for d in [intervals_dir, read_counts_dir, ploidy_dir, model_dir, calls_dir, results_dir]:
            os.makedirs(d, exist_ok=True)
        
        # Initialize Docker runner
        docker_runner = DockerRunner(device=device)
        
        # Get workdir for Docker (parent of output_dir)
        workdir = str(Path(output_dir).absolute().parent.parent)
        
        # Step 1: PreprocessIntervals
        print("Step 1/7: PreprocessIntervals...")
        preprocessed_intervals = os.path.join(intervals_dir, "preprocessed.interval_list")
        caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS] = preprocessed_intervals
        
        if os.path.exists(preprocessed_intervals):
            print(f"  Output already exists, skipping: {preprocessed_intervals}")
        else:
            preprocess_cmd = cls._create_preprocess_intervals_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK,
                " ".join(preprocess_cmd),
                workdir=workdir,
            )
            if return_code not in (0, None):
                if os.path.exists(preprocessed_intervals):
                    os.remove(preprocessed_intervals)
                raise RuntimeError(f"PreprocessIntervals failed with return code {return_code}")
        
        # Step 2: AnnotateIntervals
        print("Step 2/7: AnnotateIntervals...")
        annotated_intervals = os.path.join(intervals_dir, "annotated.tsv")
        caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS] = annotated_intervals
        
        if os.path.exists(annotated_intervals):
            print(f"  Output already exists, skipping: {annotated_intervals}")
        else:
            annotate_cmd = cls._create_annotate_intervals_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK,
                " ".join(annotate_cmd),
                workdir=workdir,
            )
            if return_code not in (0, None):
                if os.path.exists(annotated_intervals):
                    os.remove(annotated_intervals)
                raise RuntimeError(f"AnnotateIntervals failed with return code {return_code}")
        
        # Step 3: CollectReadCounts for each sample
        print("Step 3/7: CollectReadCounts...")
        sample_inputs = caller_config[GermlineCNVCallingKeys.SAMPLE_INPUTS]
        read_counts = {}
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        
        for i, bam_file in enumerate(sample_inputs):
            sample_name = os.path.basename(bam_file).replace(".bam", "")
            rc_output = os.path.join(read_counts_dir, f"read_counts_{sample_name}.hdf5")
            read_counts[sample_name] = rc_output
            
            if os.path.exists(rc_output):
                print(f"  Output already exists, skipping: {sample_name}")
            else:
                print(f"  Collecting read counts for {sample_name}...")
                collect_cmd = cls._create_collect_read_counts_command(
                    bam_file,
                    preprocessed_intervals,
                    rc_output,
                    library_paths,
                    params,
                )
                return_code = docker_runner.run(
                    DockerImages.GATK,
                    " ".join(collect_cmd),
                    workdir=workdir,
                )
                if return_code not in (0, None):
                    if os.path.exists(rc_output):
                        os.remove(rc_output)
                    raise RuntimeError(f"CollectReadCounts failed for {sample_name} with return code {return_code}")
        
        caller_config[GermlineCNVCallingKeys.READ_COUNTS] = read_counts
        
        # Step 4: FilterIntervals
        print("Step 4/7: FilterIntervals...")
        filtered_intervals = os.path.join(intervals_dir, "filtered.interval_list")
        caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS] = filtered_intervals
        
        if os.path.exists(filtered_intervals):
            print(f"  Output already exists, skipping: {filtered_intervals}")
        else:
            filter_cmd = cls._create_filter_intervals_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK,
                " ".join(filter_cmd),
                workdir=workdir,
            )
            if return_code not in (0, None):
                if os.path.exists(filtered_intervals):
                    os.remove(filtered_intervals)
                raise RuntimeError(f"FilterIntervals failed with return code {return_code}")
        
        # Step 5: DetermineGermlineContigPloidy
        print("Step 5/7: DetermineGermlineContigPloidy...")
        caller_config[GermlineCNVCallingKeys.OUTPUT_DIR] = ploidy_dir
        
        ploidy_model_dir = os.path.join(ploidy_dir, "ploidy-model")
        ploidy_calls_dir = os.path.join(ploidy_dir, "ploidy-calls")
        
        if os.path.exists(ploidy_model_dir) and os.path.exists(ploidy_calls_dir):
            print(f"  Output already exists, skipping: {ploidy_model_dir}")
        else:
            ploidy_cmd = cls._create_determine_ploidy_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK,
                " ".join(ploidy_cmd),
                workdir=workdir,
                run_as_root=True,
            )
            if return_code not in (0, None):
                import shutil
                if os.path.exists(ploidy_model_dir):
                    shutil.rmtree(ploidy_model_dir)
                if os.path.exists(ploidy_calls_dir):
                    shutil.rmtree(ploidy_calls_dir)
                raise RuntimeError(f"DetermineGermlineContigPloidy failed with return code {return_code}")
        
        caller_config[GermlineCNVCallingKeys.PLOIDY_MODEL] = ploidy_model_dir
        caller_config[GermlineCNVCallingKeys.PLOIDY_CALLS] = ploidy_calls_dir
        
        # Step 6: GermlineCNVCaller
        print("Step 6/7: GermlineCNVCaller...")
        caller_config[GermlineCNVCallingKeys.OUTPUT_DIR] = model_dir
        
        gcnv_model_dir = os.path.join(model_dir, "gcnv-model")
        gcnv_calls_dir = os.path.join(model_dir, "gcnv-calls")
        
        if os.path.exists(gcnv_model_dir) and os.path.exists(gcnv_calls_dir):
            print(f"  Output already exists, skipping: {gcnv_model_dir}")
        else:
            gcnv_cmd = cls._create_germline_cnv_caller_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK,
                " ".join(gcnv_cmd),
                workdir=workdir,
                run_as_root=True,
            )
            if return_code not in (0, None):
                import shutil
                if os.path.exists(gcnv_model_dir):
                    shutil.rmtree(gcnv_model_dir)
                if os.path.exists(gcnv_calls_dir):
                    shutil.rmtree(gcnv_calls_dir)
                raise RuntimeError(f"GermlineCNVCaller failed with return code {return_code}")
        
        caller_config[GermlineCNVCallingKeys.CNV_MODEL] = gcnv_model_dir
        caller_config[GermlineCNVCallingKeys.CNV_CALLS] = gcnv_calls_dir
        
        # Step 7: PostprocessGermlineCNVCalls for each sample
        print("Step 7/7: PostprocessGermlineCNVCalls...")
        caller_config[GermlineCNVCallingKeys.OUTPUT_DIR] = results_dir
        
        sample_names = list(read_counts.keys())
        for i, sample_name in enumerate(sample_names):
            intervals_output = os.path.join(
                results_dir, f"genotyped_intervals_{sample_name}.vcf.gz"
            )
            segments_output = os.path.join(
                results_dir, f"genotyped_segments_{sample_name}.vcf.gz"
            )
            denoised_output = os.path.join(
                results_dir, f"denoised_copy_ratios_{sample_name}.tsv"
            )
            
            if (os.path.exists(intervals_output) and 
                os.path.exists(segments_output) and 
                os.path.exists(denoised_output)):
                print(f"  Output already exists, skipping: {sample_name}")
            else:
                print(f"  Postprocessing {sample_name}...")
                caller_config[GermlineCNVCallingKeys.GENOTYPED_INTERVALS] = intervals_output
                caller_config[GermlineCNVCallingKeys.GENOTYPED_SEGMENTS] = segments_output
                caller_config[GermlineCNVCallingKeys.DENOISED_COPY_RATIOS] = denoised_output
                
                postprocess_cmd = cls._create_postprocess_calls_command(
                    caller_config, library_paths, app_config,
                    sample_index=i, sample_name=sample_name
                )
                return_code = docker_runner.run(
                    DockerImages.GATK,
                    " ".join(postprocess_cmd),
                    workdir=workdir,
                    run_as_root=True,
                )
                if return_code not in (0, None):
                    if os.path.exists(intervals_output):
                        os.remove(intervals_output)
                    if os.path.exists(segments_output):
                        os.remove(segments_output)
                    if os.path.exists(denoised_output):
                        os.remove(denoised_output)
                    raise RuntimeError(f"PostprocessGermlineCNVCalls failed for {sample_name} with return code {return_code}")
        
        print("GermlineCNVCaller pipeline completed successfully!")
