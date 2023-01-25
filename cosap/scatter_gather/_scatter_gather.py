from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor

from .._config import AppConfig
from .._pipeline_config import PipelineBaseKeys, VariantCallingKeys
from .utils import split_bam_by_intervals, get_region_file_list
from ..pipeline_builder import VariantCaller
from itertools import repeat


class ScatterGather:
    @staticmethod
    def split_variantcaller_configs(
        config: dict, split_bams: bool = False
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
                split_bam_by_intervals(germline_bam) if germline_bam else repeat(None)
            )
            splitted_tumor_bams = (
                split_bam_by_intervals(tumor_bam) if tumor_bam else repeat(None)
            )
            bam_pairs = list(zip(splitted_germline_bams, splitted_tumor_bams))


        interval_files = get_region_file_list(file_type="interval_list")
        splitted_configs = []

        for i in range(len(interval_files)):
            unfiltered_output_file = config[
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
            ]
            snp_output_file = config[VariantCallingKeys.SNP_OUTPUT]
            indel_output_file = config[VariantCallingKeys.INDEL_OUTPUT]
            gvcf_output_file = config[VariantCallingKeys.GVCF_OUTPUT]
            other_variants_output_file = config[
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT
            ]
            all_variants_output_file = config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
            output_dir = config[VariantCallingKeys.OUTPUT_DIR]
            cnf = {
                VariantCallingKeys.GERMLINE_INPUT: bam_pairs[i][0] if split_bams else germline_bam,
                VariantCallingKeys.TUMOR_INPUT: bam_pairs[i][1] if split_bams else tumor_bam,
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: f"{unfiltered_output_file}.tmp{i}",
                VariantCallingKeys.ALL_VARIANTS_OUTPUT: f"{all_variants_output_file}.tmp{i}",
                VariantCallingKeys.SNP_OUTPUT: f"{snp_output_file}.tmp{i}",
                VariantCallingKeys.INDEL_OUTPUT: f"{indel_output_file}.tmp{i}",
                VariantCallingKeys.GVCF_OUTPUT: f"{gvcf_output_file}.tmp{i}",
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: f"{other_variants_output_file}.tmp{i}",
                VariantCallingKeys.OUTPUT_DIR: output_dir,
                VariantCallingKeys.BED_FILE: interval_files[i],
            }
            splitted_configs.append(cnf)

        return splitted_configs

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
    def run_splitted_configs(run_function: Callable, func_params: list):
        app_config = AppConfig()
        with ProcessPoolExecutor(
            max_workers=app_config.MAX_THREADS_PER_JOB
        ) as executor:
            executor.map(run_function, func_params)
