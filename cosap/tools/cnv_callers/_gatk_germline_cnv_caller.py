"""
GATK GermlineCNVCaller Pipeline Implementation

This module implements the complete GATK GermlineCNVCaller workflow for detecting
rare germline copy number variants (CNVs) as described in:
https://gatk.broadinstitute.org/hc/en-us/articles/360035531152

The workflow consists of the following steps:
1. PreprocessIntervals - Prepare intervals for analysis (bin WGS or pad exomes)
2. AnnotateIntervals - Annotate intervals with GC content and mappability (optional but recommended)
3. CollectReadCounts - Collect read counts per interval per sample (PARALLELIZED)
4. FilterIntervals - Filter intervals based on annotations and extreme counts
5. ScatterIntervals - Split intervals into shards for parallel processing (NEW)
6. DetermineGermlineContigPloidy - Determine baseline contig ploidies
7. GermlineCNVCaller - Call germline CNVs (cohort or case mode) (PARALLELIZED via sharding)
8. PostprocessGermlineCNVCalls - Generate per-sample VCF results (PARALLELIZED)

Parallelization Strategy:
- Step 3 (CollectReadCounts): Parallel across samples
- Step 7 (GermlineCNVCaller): Parallel across interval shards
- Step 8 (PostprocessGermlineCNVCalls): Parallel across samples
"""

import glob
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import GermlineCNVCallingKeys
from ...runners.docker_runner import DockerRunner
from ._cnv_callers import _CNVCaller


@dataclass
class ParallelTaskResult:
    """Result of a parallel task execution."""

    task_id: str
    output_path: str
    success: bool
    error_message: Optional[str] = None


class GATKGermlineCNVCaller(_CNVCaller):
    """
    GATK GermlineCNVCaller - Complete pipeline for germline CNV calling.

    This class orchestrates the entire GATK GermlineCNVCaller workflow,
    executing all steps in sequence to produce final VCF outputs.
    All GATK commands are executed inside Docker containers.

    Parallelization:
    - CollectReadCounts runs in parallel across samples
    - GermlineCNVCaller runs in parallel across interval shards
    - PostprocessGermlineCNVCalls runs in parallel across samples
    """

    # Default intervals per shard (~5Mbp for 1kb bins)
    DEFAULT_INTERVALS_PER_SHARD = 5000

    @classmethod
    def _get_standard_contigs(cls, params: Dict) -> List[str]:
        """Get list of standard contigs (chr1-22, chrX, chrY) based on reference style."""
        include_contigs = params.get("include_contigs", None)
        if include_contigs:
            return include_contigs

        use_chr_prefix = params.get("use_chr_prefix", True)
        if use_chr_prefix:
            return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        else:
            return [str(i) for i in range(1, 23)] + ["X", "Y"]

    @classmethod
    def _create_preprocess_intervals_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """Create command for GATK PreprocessIntervals."""
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        bin_length = params.get("bin_length", 1000)
        padding = params.get("padding", 0)

        command = [
            "gatk",
            "PreprocessIntervals",
            "-R",
            library_paths.REF_FASTA,
            "--bin-length",
            str(bin_length),
            "--padding",
            str(padding),
            "--interval-merging-rule",
            "OVERLAPPING_ONLY",
            "-O",
            caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS],
        ]

        if GermlineCNVCallingKeys.INTERVALS in caller_config:
            command.extend(["-L", caller_config[GermlineCNVCallingKeys.INTERVALS]])

        standard_contigs = cls._get_standard_contigs(params)
        for contig in standard_contigs:
            command.extend(["-L", contig])

        if params.get("exclude_intervals"):
            exclude = params["exclude_intervals"]
            if isinstance(exclude, list):
                for interval in exclude:
                    command.extend(["-XL", interval])
            else:
                command.extend(["-XL", exclude])

        return command

    @classmethod
    def _create_annotate_intervals_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """Create command for GATK AnnotateIntervals."""
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})

        command = [
            "gatk",
            "AnnotateIntervals",
            "-R",
            library_paths.REF_FASTA,
            "-L",
            caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS],
            "--interval-merging-rule",
            "OVERLAPPING_ONLY",
            "-O",
            caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS],
        ]

        if params.get("mappability_track"):
            command.extend(["--mappability-track", params["mappability_track"]])
        if params.get("segmental_duplication_track"):
            command.extend(
                ["--segmental-duplication-track", params["segmental_duplication_track"]]
            )

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
        """Create command for GATK CollectReadCounts."""
        params = params or {}

        command = [
            "gatk",
            "CollectReadCounts",
            "-I",
            bam_file,
            "-L",
            intervals,
            "-R",
            library_paths.REF_FASTA,
            "--interval-merging-rule",
            "OVERLAPPING_ONLY",
            "-O",
            output_file,
        ]

        if params.get("format") == "TSV":
            command.extend(["--format", "TSV"])

        return command

    @classmethod
    def _create_scatter_intervals_command(
        cls,
        input_intervals: str,
        output_dir: str,
        intervals_per_shard: int,
    ) -> List[str]:
        """
        Create command for Picard IntervalListTools to scatter intervals.

        This splits the filtered intervals into shards for parallel GermlineCNVCaller execution.
        Each shard should have ~5000 intervals (~5Mbp for 1kb bins) as recommended by GATK.
        """
        command = [
            "gatk",
            "IntervalListTools",
            "--INPUT",
            input_intervals,
            "--SUBDIVISION_MODE",
            "INTERVAL_COUNT",
            "--SCATTER_CONTENT",
            str(intervals_per_shard),
            "--OUTPUT",
            output_dir,
        ]
        return command

    @classmethod
    def _create_filter_intervals_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List[str]:
        """Create command for GATK FilterIntervals."""
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})

        command = [
            "gatk",
            "FilterIntervals",
            "-L",
            caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS],
            "--annotated-intervals",
            caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS],
            "--interval-merging-rule",
            "OVERLAPPING_ONLY",
            "-O",
            caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS],
        ]

        read_counts = caller_config[GermlineCNVCallingKeys.READ_COUNTS]
        if isinstance(read_counts, dict):
            for rc_file in read_counts.values():
                command.extend(["-I", rc_file])
        else:
            for rc_file in read_counts:
                command.extend(["-I", rc_file])

        # Filter thresholds
        threshold_params = [
            ("minimum_gc_content", "--minimum-gc-content"),
            ("maximum_gc_content", "--maximum-gc-content"),
            ("minimum_mappability", "--minimum-mappability"),
            ("maximum_mappability", "--maximum-mappability"),
            ("low_count_filter_count_threshold", "--low-count-filter-count-threshold"),
            (
                "low_count_filter_percentage_of_samples",
                "--low-count-filter-percentage-of-samples",
            ),
            (
                "extreme_count_filter_minimum_percentile",
                "--extreme-count-filter-minimum-percentile",
            ),
            (
                "extreme_count_filter_maximum_percentile",
                "--extreme-count-filter-maximum-percentile",
            ),
            (
                "extreme_count_filter_percentage_of_samples",
                "--extreme-count-filter-percentage-of-samples",
            ),
        ]

        for param_key, flag in threshold_params:
            if params.get(param_key):
                command.extend([flag, str(params[param_key])])

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
        """Create command for GATK DetermineGermlineContigPloidy."""
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        run_mode = caller_config.get(GermlineCNVCallingKeys.RUN_MODE, "COHORT")

        command = [
            "gatk",
            "DetermineGermlineContigPloidy",
            "-L",
            caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS],
            "--interval-merging-rule",
            "OVERLAPPING_ONLY",
            "--output",
            caller_config[GermlineCNVCallingKeys.OUTPUT_DIR],
            "--output-prefix",
            "ploidy",
        ]

        read_counts = caller_config[GermlineCNVCallingKeys.READ_COUNTS]
        if isinstance(read_counts, dict):
            for rc_file in read_counts.values():
                command.extend(["-I", rc_file])
        else:
            for rc_file in read_counts:
                command.extend(["-I", rc_file])

        if run_mode == "COHORT":
            if GermlineCNVCallingKeys.CONTIG_PLOIDY_PRIORS in caller_config:
                command.extend(
                    [
                        "--contig-ploidy-priors",
                        caller_config[GermlineCNVCallingKeys.CONTIG_PLOIDY_PRIORS],
                    ]
                )
            elif hasattr(library_paths, "CONTIG_PLOIDY_PRIORS"):
                command.extend(
                    ["--contig-ploidy-priors", library_paths.CONTIG_PLOIDY_PRIORS]
                )
        else:
            if params.get("ploidy_model"):
                command.extend(["--model", params["ploidy_model"]])

        command.extend(["--verbosity", params.get("verbosity", "INFO")])
        return command

    @classmethod
    def _create_germline_cnv_caller_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        shard_intervals: str,
        output_prefix: str,
    ) -> List[str]:
        """
        Create command for GATK GermlineCNVCaller for a single shard.

        Each shard processes a subset of intervals independently.
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})
        run_mode = caller_config.get(GermlineCNVCallingKeys.RUN_MODE, "COHORT")

        command = [
            "gatk",
            "GermlineCNVCaller",
            "--run-mode",
            run_mode,
            "-L",
            shard_intervals,
            "--contig-ploidy-calls",
            caller_config[GermlineCNVCallingKeys.PLOIDY_CALLS],
            "--interval-merging-rule",
            "OVERLAPPING_ONLY",
            "--output",
            caller_config[GermlineCNVCallingKeys.OUTPUT_DIR],
            "--output-prefix",
            output_prefix,
        ]

        read_counts = caller_config[GermlineCNVCallingKeys.READ_COUNTS]
        if isinstance(read_counts, dict):
            for rc_file in read_counts.values():
                command.extend(["-I", rc_file])
        else:
            for rc_file in read_counts:
                command.extend(["-I", rc_file])

        if GermlineCNVCallingKeys.ANNOTATED_INTERVALS in caller_config:
            command.extend(
                [
                    "--annotated-intervals",
                    caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS],
                ]
            )

        if run_mode == "CASE" and params.get("cnv_model"):
            command.extend(["--model", params["cnv_model"]])

        # Sensitivity tuning parameters
        sensitivity_params = [
            ("class_coherence_length", "--class-coherence-length"),
            ("cnv_coherence_length", "--cnv-coherence-length"),
            ("interval_psi_scale", "--interval-psi-scale"),
            ("log_mean_bias_standard_deviation", "--log-mean-bias-standard-deviation"),
            ("sample_psi_scale", "--sample-psi-scale"),
            ("depth_correction_tau", "--depth-correction-tau"),
            ("p_active", "--p-active"),
            ("p_alt", "--p-alt"),
            ("max_copy_number", "--max-copy-number"),
        ]

        for param_key, flag in sensitivity_params:
            if params.get(param_key):
                command.extend([flag, str(params[param_key])])

        if params.get("enable_bias_factors") is not None:
            command.extend(
                ["--enable-bias-factors", str(params["enable_bias_factors"]).lower()]
            )

        command.extend(["--verbosity", params.get("verbosity", "INFO")])
        return command

    @classmethod
    def _create_postprocess_calls_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        sample_index: int,
        intervals_output: str,
        segments_output: str,
        denoised_output: str,
        model_shards: List[str],
        calls_shards: List[str],
    ) -> List[str]:
        """
        Create command for GATK PostprocessGermlineCNVCalls.

        Consolidates results from all shards for a single sample.
        """
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})

        if hasattr(library_paths, "REF_DICT"):
            sequence_dict = library_paths.REF_DICT
        else:
            sequence_dict = library_paths.REF_FASTA.replace(".fasta", ".dict").replace(
                ".fa", ".dict"
            )

        command = [
            "gatk",
            "PostprocessGermlineCNVCalls",
            "--sample-index",
            str(sample_index),
            "--output-genotyped-intervals",
            intervals_output,
            "--output-genotyped-segments",
            segments_output,
            "--output-denoised-copy-ratios",
            denoised_output,
            "--sequence-dictionary",
            sequence_dict,
        ]

        # Add all model shards (order matters!)
        for shard in sorted(model_shards):
            command.extend(["--model-shard-path", shard])

        # Add all calls shards (order matters!)
        for shard in sorted(calls_shards):
            command.extend(["--calls-shard-path", shard])

        command.extend(
            [
                "--contig-ploidy-calls",
                caller_config[GermlineCNVCallingKeys.PLOIDY_CALLS],
            ]
        )

        allosomal_contigs = params.get("allosomal_contigs", ["chrX", "chrY"])
        for contig in allosomal_contigs:
            command.extend(["--allosomal-contig", contig])

        if params.get("autosomal_ref_copy_number"):
            command.extend(
                [
                    "--autosomal-ref-copy-number",
                    str(params["autosomal_ref_copy_number"]),
                ]
            )

        return command

    # ==================== PARALLEL EXECUTION METHODS ====================

    @classmethod
    def _run_docker_task(
        cls,
        command: List[str],
        task_id: str,
        output_path: str,
        workdir: str,
        device: str,
        run_as_root: bool = False,
    ) -> ParallelTaskResult:
        """Execute a single Docker task. Thread-safe by creating its own DockerRunner."""
        try:
            docker_runner = DockerRunner(device=device)
            return_code = docker_runner.run(
                DockerImages.GATK,
                " ".join(command),
                workdir=workdir,
                run_as_root=run_as_root,
            )

            if return_code not in (0, None):
                return ParallelTaskResult(
                    task_id=task_id,
                    output_path=output_path,
                    success=False,
                    error_message=f"Command failed with return code {return_code}",
                )

            return ParallelTaskResult(
                task_id=task_id,
                output_path=output_path,
                success=True,
            )
        except Exception as e:
            return ParallelTaskResult(
                task_id=task_id,
                output_path=output_path,
                success=False,
                error_message=str(e),
            )

    @classmethod
    def _collect_read_counts_parallel(
        cls,
        sample_inputs: List[str],
        preprocessed_intervals: str,
        read_counts_dir: str,
        library_paths: LibraryPaths,
        params: Dict,
        workdir: str,
        device: str,
        max_workers: Optional[int] = None,
    ) -> Dict[str, str]:
        """Collect read counts for all samples in parallel."""
        tasks = []
        for bam_file in sample_inputs:
            sample_name = os.path.basename(bam_file).replace(".bam", "")
            rc_output = os.path.join(read_counts_dir, f"read_counts_{sample_name}.hdf5")

            if os.path.exists(rc_output):
                tasks.append((sample_name, rc_output, None))  # Already exists
            else:
                command = cls._create_collect_read_counts_command(
                    bam_file, preprocessed_intervals, rc_output, library_paths, params
                )
                tasks.append((sample_name, rc_output, command))

        if max_workers is None:
            max_workers = min(len(tasks), os.cpu_count() or 4)

        read_counts = {}
        failed = []

        # Separate existing from needing processing
        to_process = [
            (name, output, cmd) for name, output, cmd in tasks if cmd is not None
        ]
        for name, output, _ in tasks:
            if os.path.exists(output):
                read_counts[name] = output
                print(f"    ✓ {name} (cached)")

        if to_process:
            print(
                f"  Processing {len(to_process)} samples with {max_workers} workers..."
            )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        cls._run_docker_task, cmd, name, output, workdir, device
                    ): name
                    for name, output, cmd in to_process
                }

                for future in as_completed(futures):
                    result = future.result()
                    if result.success:
                        read_counts[result.task_id] = result.output_path
                        print(f"    ✓ {result.task_id}")
                    else:
                        failed.append((result.task_id, result.error_message))
                        print(f"    ✗ {result.task_id}: {result.error_message}")
                        if os.path.exists(result.output_path):
                            os.remove(result.output_path)

        if failed:
            raise RuntimeError(
                f"CollectReadCounts failed for: {[f[0] for f in failed]}"
            )

        return read_counts

    @classmethod
    def _scatter_intervals(
        cls,
        filtered_intervals: str,
        scatter_dir: str,
        intervals_per_shard: int,
        workdir: str,
        device: str,
    ) -> List[str]:
        """
        Scatter intervals into shards using IntervalListTools.

        Returns list of paths to scattered interval files.
        """
        os.makedirs(scatter_dir, exist_ok=True)

        # Check if already scattered
        existing_shards = glob.glob(
            os.path.join(scatter_dir, "*/scattered.interval_list")
        )
        if existing_shards:
            print(f"  Using {len(existing_shards)} existing shards")
            return sorted(existing_shards)

        command = cls._create_scatter_intervals_command(
            filtered_intervals, scatter_dir, intervals_per_shard
        )

        docker_runner = DockerRunner(device=device)
        return_code = docker_runner.run(
            DockerImages.GATK,
            " ".join(command),
            workdir=workdir,
        )

        if return_code not in (0, None):
            raise RuntimeError(
                f"IntervalListTools scatter failed with return code {return_code}"
            )

        # Find scattered intervals (IntervalListTools creates temp_XXXX_of_N directories)
        shards = sorted(
            glob.glob(os.path.join(scatter_dir, "*/scattered.interval_list"))
        )

        if not shards:
            raise RuntimeError(f"No scattered intervals found in {scatter_dir}")

        print(f"  Created {len(shards)} interval shards")
        return shards

    @classmethod
    def _run_gcnv_caller_parallel(
        cls,
        shard_intervals: List[str],
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        model_dir: str,
        workdir: str,
        device: str,
        max_workers: Optional[int] = None,
    ) -> Tuple[List[str], List[str]]:
        """
        Run GermlineCNVCaller in parallel across interval shards.

        Returns:
            Tuple of (model_shard_paths, calls_shard_paths)
        """
        tasks = []

        for i, shard_interval in enumerate(shard_intervals):
            shard_name = f"shard_{i:04d}"
            output_prefix = shard_name

            model_output = os.path.join(model_dir, f"{output_prefix}-model")
            calls_output = os.path.join(model_dir, f"{output_prefix}-calls")

            if os.path.exists(model_output) and os.path.exists(calls_output):
                tasks.append((shard_name, model_output, calls_output, None))
            else:
                shard_config = caller_config.copy()
                shard_config[GermlineCNVCallingKeys.OUTPUT_DIR] = model_dir

                command = cls._create_germline_cnv_caller_command(
                    shard_config,
                    library_paths,
                    app_config,
                    shard_interval,
                    output_prefix,
                )
                tasks.append((shard_name, model_output, calls_output, command))

        if max_workers is None:
            max_workers = min(len(tasks), os.cpu_count() or 4)

        model_shards = []
        calls_shards = []
        failed = []

        # Collect cached results
        to_process = []
        for shard_name, model_out, calls_out, cmd in tasks:
            if cmd is None:
                model_shards.append(model_out)
                calls_shards.append(calls_out)
                print(f"    ✓ {shard_name} (cached)")
            else:
                to_process.append((shard_name, model_out, calls_out, cmd))

        if to_process:
            print(
                f"  Processing {len(to_process)} shards with {max_workers} workers..."
            )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        cls._run_docker_task,
                        cmd,
                        shard_name,
                        model_out,
                        workdir,
                        device,
                        True,
                    ): (shard_name, model_out, calls_out)
                    for shard_name, model_out, calls_out, cmd in to_process
                }

                for future in as_completed(futures):
                    shard_name, model_out, calls_out = futures[future]
                    result = future.result()

                    if result.success:
                        model_shards.append(model_out)
                        calls_shards.append(calls_out)
                        print(f"    ✓ {shard_name}")
                    else:
                        failed.append((shard_name, result.error_message))
                        print(f"    ✗ {shard_name}: {result.error_message}")
                        # Cleanup partial outputs
                        if os.path.exists(model_out):
                            shutil.rmtree(model_out)
                        if os.path.exists(calls_out):
                            shutil.rmtree(calls_out)

        if failed:
            raise RuntimeError(
                f"GermlineCNVCaller failed for shards: {[f[0] for f in failed]}"
            )

        return sorted(model_shards), sorted(calls_shards)

    @classmethod
    def _postprocess_calls_parallel(
        cls,
        sample_names: List[str],
        model_shards: List[str],
        calls_shards: List[str],
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        results_dir: str,
        workdir: str,
        device: str,
        max_workers: Optional[int] = None,
    ) -> List[Tuple[str, str, str]]:
        """Postprocess CNV calls for all samples in parallel."""
        tasks = []

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

            if all(
                os.path.exists(f)
                for f in [intervals_output, segments_output, denoised_output]
            ):
                tasks.append(
                    (
                        sample_name,
                        i,
                        intervals_output,
                        segments_output,
                        denoised_output,
                        None,
                    )
                )
            else:
                command = cls._create_postprocess_calls_command(
                    caller_config,
                    library_paths,
                    app_config,
                    sample_index=i,
                    intervals_output=intervals_output,
                    segments_output=segments_output,
                    denoised_output=denoised_output,
                    model_shards=model_shards,
                    calls_shards=calls_shards,
                )
                tasks.append(
                    (
                        sample_name,
                        i,
                        intervals_output,
                        segments_output,
                        denoised_output,
                        command,
                    )
                )

        if max_workers is None:
            max_workers = min(len(tasks), os.cpu_count() or 4)

        results = []
        failed = []

        to_process = []
        for sample_name, idx, int_out, seg_out, den_out, cmd in tasks:
            if cmd is None:
                results.append((sample_name, int_out, seg_out))
                print(f"    ✓ {sample_name} (cached)")
            else:
                to_process.append((sample_name, idx, int_out, seg_out, den_out, cmd))

        if to_process:
            print(
                f"  Processing {len(to_process)} samples with {max_workers} workers..."
            )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        cls._run_docker_task,
                        cmd,
                        sample_name,
                        seg_out,
                        workdir,
                        device,
                        True,
                    ): (sample_name, int_out, seg_out, den_out)
                    for sample_name, _, int_out, seg_out, den_out, cmd in to_process
                }

                for future in as_completed(futures):
                    sample_name, int_out, seg_out, den_out = futures[future]
                    result = future.result()

                    if result.success:
                        results.append((sample_name, int_out, seg_out))
                        print(f"    ✓ {sample_name}")
                    else:
                        failed.append((sample_name, result.error_message))
                        print(f"    ✗ {sample_name}: {result.error_message}")
                        for f in [int_out, seg_out, den_out]:
                            if os.path.exists(f):
                                os.remove(f)

        if failed:
            raise RuntimeError(
                f"PostprocessGermlineCNVCalls failed for: {[f[0] for f in failed]}"
            )

        return results

    # ==================== MAIN PIPELINE ====================

    @classmethod
    def call(
        cls, caller_config: Dict, device: str = "cpu", max_workers: Optional[int] = None
    ):
        """
        Execute the complete GATK GermlineCNVCaller pipeline with parallelization.

        Pipeline steps:
        1. PreprocessIntervals - Prepare analysis intervals
        2. AnnotateIntervals - Add GC content annotations
        3. CollectReadCounts - Collect per-sample read counts (PARALLEL)
        4. FilterIntervals - Filter problematic intervals
        5. ScatterIntervals - Split into shards for parallel processing
        6. DetermineGermlineContigPloidy - Determine contig ploidies
        7. GermlineCNVCaller - Call CNVs (PARALLEL across shards)
        8. PostprocessGermlineCNVCalls - Generate final VCFs (PARALLEL)

        Args:
            caller_config: Configuration dictionary
            device: "cpu" or "gpu"
            max_workers: Max parallel workers (default: CPU count)
        """
        library_paths = LibraryPaths()
        app_config = AppConfig()
        params = caller_config.get(GermlineCNVCallingKeys.PARAMS, {})

        output_dir = caller_config[GermlineCNVCallingKeys.OUTPUT_DIR]
        os.makedirs(output_dir, exist_ok=True)

        # Create subdirectories
        intervals_dir = os.path.join(output_dir, "intervals")
        scatter_dir = os.path.join(output_dir, "scatter")
        read_counts_dir = os.path.join(output_dir, "read_counts")
        ploidy_dir = os.path.join(output_dir, "ploidy")
        model_dir = os.path.join(output_dir, "model")
        results_dir = os.path.join(output_dir, "results")

        for d in [
            intervals_dir,
            scatter_dir,
            read_counts_dir,
            ploidy_dir,
            model_dir,
            results_dir,
        ]:
            os.makedirs(d, exist_ok=True)

        docker_runner = DockerRunner(device=device)
        workdir = str(Path(caller_config[GermlineCNVCallingKeys.WORKDIR]).absolute())

        # Step 1: PreprocessIntervals
        print("Step 1/8: PreprocessIntervals...")
        preprocessed_intervals = os.path.join(
            intervals_dir, "preprocessed.interval_list"
        )
        caller_config[GermlineCNVCallingKeys.PREPROCESSED_INTERVALS] = (
            preprocessed_intervals
        )

        if not os.path.exists(preprocessed_intervals):
            cmd = cls._create_preprocess_intervals_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK, " ".join(cmd), workdir=workdir
            )
            if return_code not in (0, None):
                raise RuntimeError(f"PreprocessIntervals failed: {return_code}")
        else:
            print(f"  Cached: {preprocessed_intervals}")

        # Step 2: AnnotateIntervals
        print("Step 2/8: AnnotateIntervals...")
        annotated_intervals = os.path.join(intervals_dir, "annotated.tsv")
        caller_config[GermlineCNVCallingKeys.ANNOTATED_INTERVALS] = annotated_intervals

        if not os.path.exists(annotated_intervals):
            cmd = cls._create_annotate_intervals_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK, " ".join(cmd), workdir=workdir
            )
            if return_code not in (0, None):
                raise RuntimeError(f"AnnotateIntervals failed: {return_code}")
        else:
            print(f"  Cached: {annotated_intervals}")

        # Step 3: CollectReadCounts (PARALLEL)
        print("Step 3/8: CollectReadCounts (parallel)...")
        sample_inputs = caller_config[GermlineCNVCallingKeys.SAMPLE_INPUTS]
        read_counts = cls._collect_read_counts_parallel(
            sample_inputs,
            preprocessed_intervals,
            read_counts_dir,
            library_paths,
            params,
            workdir,
            device,
            max_workers,
        )
        caller_config[GermlineCNVCallingKeys.READ_COUNTS] = read_counts

        # Step 4: FilterIntervals
        print("Step 4/8: FilterIntervals...")
        filtered_intervals = os.path.join(intervals_dir, "filtered.interval_list")
        caller_config[GermlineCNVCallingKeys.FILTERED_INTERVALS] = filtered_intervals

        if not os.path.exists(filtered_intervals):
            cmd = cls._create_filter_intervals_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK, " ".join(cmd), workdir=workdir
            )
            if return_code not in (0, None):
                raise RuntimeError(f"FilterIntervals failed: {return_code}")
        else:
            print(f"  Cached: {filtered_intervals}")

        # Step 5: ScatterIntervals
        print("Step 5/8: ScatterIntervals...")
        intervals_per_shard = params.get(
            "intervals_per_shard", cls.DEFAULT_INTERVALS_PER_SHARD
        )
        shard_intervals = cls._scatter_intervals(
            filtered_intervals, scatter_dir, intervals_per_shard, workdir, device
        )

        # Step 6: DetermineGermlineContigPloidy
        print("Step 6/8: DetermineGermlineContigPloidy...")
        caller_config[GermlineCNVCallingKeys.OUTPUT_DIR] = ploidy_dir

        ploidy_model_dir = os.path.join(ploidy_dir, "ploidy-model")
        ploidy_calls_dir = os.path.join(ploidy_dir, "ploidy-calls")

        if not (os.path.exists(ploidy_model_dir) and os.path.exists(ploidy_calls_dir)):
            cmd = cls._create_determine_ploidy_command(
                caller_config, library_paths, app_config
            )
            return_code = docker_runner.run(
                DockerImages.GATK, " ".join(cmd), workdir=workdir, run_as_root=True
            )
            if return_code not in (0, None):
                raise RuntimeError(
                    f"DetermineGermlineContigPloidy failed: {return_code}"
                )
        else:
            print(f"  Cached: {ploidy_model_dir}")

        caller_config[GermlineCNVCallingKeys.PLOIDY_MODEL] = ploidy_model_dir
        caller_config[GermlineCNVCallingKeys.PLOIDY_CALLS] = ploidy_calls_dir

        # Step 7: GermlineCNVCaller (PARALLEL across shards)
        print("Step 7/8: GermlineCNVCaller (parallel shards)...")
        model_shards, calls_shards = cls._run_gcnv_caller_parallel(
            shard_intervals,
            caller_config,
            library_paths,
            app_config,
            model_dir,
            workdir,
            device,
            max_workers,
        )

        caller_config[GermlineCNVCallingKeys.CNV_MODEL] = model_shards
        caller_config[GermlineCNVCallingKeys.CNV_CALLS] = calls_shards

        # Step 8: PostprocessGermlineCNVCalls (PARALLEL)
        print("Step 8/8: PostprocessGermlineCNVCalls (parallel)...")
        sample_names = list(read_counts.keys())

        cls._postprocess_calls_parallel(
            sample_names,
            model_shards,
            calls_shards,
            caller_config,
            library_paths,
            app_config,
            results_dir,
            workdir,
            device,
            max_workers,
        )

        print("=" * 60)
        print("GermlineCNVCaller pipeline completed successfully!")
        print(f"  Samples processed: {len(sample_names)}")
        print(f"  Shards used: {len(shard_intervals)}")
        print(f"  Results in: {results_dir}")
