from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor

from .._config import AppConfig
from .._pipeline_config import PipelineBaseKeys, VariantCallingKeys
from .utils import split_bam_by_intervals, get_region_file_list, create_tmp_filename
from itertools import repeat
from subprocess import run
from itertools import chain
import os


class ScatterGather:
    @staticmethod
    def split_variantcaller_configs(
        config: dict, bed_file, split_bams: bool = False
    ) -> list[dict]:
        germline_bam = (
            config[VariantCallingKeys.GERMLINE_INPUT]
            if VariantCallingKeys.GERMLINE_INPUT in config.keys()
            else None
        )
        tumor_bam = (
            config[VariantCallingKeys.TUMOR_INPUT]
            if VariantCallingKeys.TUMOR_INPUT in config.keys()
            else None
        )

        if split_bams:
            splitted_germline_bams = (
                split_bam_by_intervals(germline_bam, bed_file=bed_file)
                if germline_bam
                else repeat(None)
            )
            splitted_tumor_bams = (
                split_bam_by_intervals(tumor_bam, bed_file=bed_file)
                if tumor_bam
                else repeat(None)
            )
            bam_pairs = list(zip(splitted_germline_bams, splitted_tumor_bams))

        interval_files = get_region_file_list(
            file_type="interval_list", bed_file=config[VariantCallingKeys.BED_FILE]
        )
        splitted_configs = []

        for i in range(len(interval_files)):
            unfiltered_output_file = create_tmp_filename(config[
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
            ],i)
            snp_output_file = create_tmp_filename(
                config[VariantCallingKeys.SNP_OUTPUT], i
            )
            indel_output_file = create_tmp_filename(
                config[VariantCallingKeys.INDEL_OUTPUT], i
            )
            gvcf_output_file = create_tmp_filename(
                config[VariantCallingKeys.GVCF_OUTPUT], i
            )
            other_variants_output_file = create_tmp_filename(
                config[VariantCallingKeys.OTHER_VARIANTS_OUTPUT], i
            )
            all_variants_output_file = create_tmp_filename(
                config[VariantCallingKeys.ALL_VARIANTS_OUTPUT], i
            )
            output_dir = config[VariantCallingKeys.OUTPUT_DIR]

            cnf = {
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: unfiltered_output_file,
                VariantCallingKeys.ALL_VARIANTS_OUTPUT: all_variants_output_file,
                VariantCallingKeys.SNP_OUTPUT: snp_output_file,
                VariantCallingKeys.INDEL_OUTPUT: indel_output_file,
                VariantCallingKeys.GVCF_OUTPUT: gvcf_output_file,
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: other_variants_output_file,
                VariantCallingKeys.OUTPUT_DIR: output_dir,
                VariantCallingKeys.BED_FILE: interval_files[i],
            }

            if VariantCallingKeys.GERMLINE_INPUT in config.keys():
                cnf[VariantCallingKeys.GERMLINE_INPUT] = (
                    bam_pairs[i][0] if split_bams else germline_bam
                )
            if VariantCallingKeys.TUMOR_INPUT in config.keys():
                cnf[VariantCallingKeys.TUMOR_INPUT] = (
                    bam_pairs[i][1] if split_bams else tumor_bam
                )
            splitted_configs.append(cnf)

        return splitted_configs

    @staticmethod
    def gather_vcfs(configs: list, output_path):
        vcfs = [cfg[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT] for cfg in configs]

        command = ["gatk", "MergeVcfs", "-O", output_path]
        command.extend(list(chain(*zip(repeat("-I"), vcfs))))
        print(command)
        run(command)

    @staticmethod
    def split_bam_process_configs(config: dict) -> list[dict]:
        splitted_configs = []

        input_bam = config[PipelineBaseKeys.INPUT]
        splitted_bams = split_bam_by_intervals(input_bam)
        for bam in splitted_bams:
            cnf = {}
            cnf[PipelineBaseKeys.INPUT] = bam
            splitted_configs.append(cnf)

        return splitted_configs

    @staticmethod
    def run_parallel(run_function: Callable, func_params: list):
        app_config = AppConfig()
        with ProcessPoolExecutor(
            max_workers=app_config.MAX_THREADS_PER_JOB
        ) as executor:
            executor.map(run_function, func_params)
