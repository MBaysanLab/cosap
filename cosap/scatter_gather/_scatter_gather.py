import glob
import os
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from itertools import chain, repeat
from subprocess import run
import shortuuid
from ..pipeline_builder import VariantCaller

from .._config import AppConfig
from .._pipeline_config import PipelineBaseKeys, VariantCallingKeys, PipelineKeys
from .utils import (create_tmp_filename, get_region_file_list,
                    split_bam_by_intervals)


class ScatterGather:
    @staticmethod
    def split_variantcaller_configs(
        config: dict, bed_file=None, split_bams: bool = False
    ) -> list[dict]:

        # If the number of threads is not suitable for parellelization, return the original config
        app_config = AppConfig()
        threads = app_config.MAX_THREADS_PER_JOB

        if int(threads) == 1:
            return [config]

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
            file_type="interval_list", bed_file=bed_file
        )
        splitted_configs = []

        for i in range(len(interval_files)):
            tmp_name = f"tmp{shortuuid.uuid()}"
            variant_caller = VariantCaller(
                library=config[VariantCallingKeys.LIBRARY],
                name=tmp_name,
                bed_file=interval_files[i],
                params=config[PipelineBaseKeys.PARAMS],
                gvcf=config[VariantCallingKeys.OUTPUT_TYPE] == "GVCF",
            )
            cfg = variant_caller.get_config()[PipelineKeys.VARIANT_CALLING][tmp_name]

            if VariantCallingKeys.GERMLINE_INPUT in config.keys():
                cfg[VariantCallingKeys.GERMLINE_INPUT] = (
                    bam_pairs[i][0] if split_bams else germline_bam
                )
            if VariantCallingKeys.TUMOR_INPUT in config.keys():
                cfg[VariantCallingKeys.TUMOR_INPUT] = (
                    bam_pairs[i][1] if split_bams else tumor_bam
                )
            cfg = dict(config, **cfg)
            splitted_configs.append(cfg)

        return splitted_configs

    @staticmethod
    def gather_vcfs(configs: list, output_path, mode="vcf"):

        GVCF_MODE = "gvcf"
        VCF_MODE = "vcf"

        if mode.lower() == VCF_MODE:
            vcfs = [cfg[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT] for cfg in configs]
            command = ["gatk", "MergeVcfs", "-O", output_path]
        elif mode.lower() == GVCF_MODE:
            vcfs = [cfg[VariantCallingKeys.GVCF_OUTPUT] for cfg in configs]
            command = ["gatk", "SortVcf", "-O", output_path]
        else:
            raise ValueError(f"Mode {mode} not supported")
        
        command.extend(list(chain(*zip(repeat("-I"), vcfs))))
        run(command)

    @staticmethod
    def clean_temp_files(path):
        temp_files = glob.glob(f"{path}/*tmp*")
        for tmp in temp_files:
            os.remove(tmp)

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
