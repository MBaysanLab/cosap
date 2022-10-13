from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._pipeline_config import PipelineBaseKeys, VariantCallingKeys
from .utils import split_bam_by_intervals


class ScatterGather:
    @staticmethod
    def split_variantcaller_configs(config: Dict) -> List[Dict]:
        splitted_configs = []

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

        if germline_bam:
            splitted_germline_bams = split_bam_by_intervals(germline_bam)
        if tumor_bam:
            splitted_tumor_bams = split_bam_by_intervals(tumor_bam)

        if germline_bam and tumor_bam:
            pairs = zip(splitted_germline_bams, splitted_tumor_bams)
            for pair in pairs:
                cnf = {}
                cnf[VariantCallingKeys.GERMLINE_INPUT] = pair[0]
                cnf[VariantCallingKeys.TUMOR_INPUT] = pair[1]
                splitted_configs.append(cnf)

        elif germline_bam:
            for bam in splitted_germline_bams:
                cnf = {}
                cnf[VariantCallingKeys.GERMLINE_INPUT] = bam
                splitted_configs.append(cnf)

        elif tumor_bam:
            for bam in splitted_tumor_bams:
                cnf = {}
                cnf[VariantCallingKeys.TUMOR_INPUT] = bam
                splitted_configs.append(cnf)

        return splitted_configs

    @staticmethod
    def split_bam_process_configs(config: Dict) -> List[Dict]:
        splitted_configs = []

        input_bam = config[PipelineBaseKeys.INPUT]
        splitted_bams = split_bam_by_intervals(input_bam)
        for bam in splitted_bams:
            cnf = {}
            cnf[PipelineBaseKeys.INPUT] = bam
            splitted_configs.append(cnf)

        return splitted_configs

    @staticmethod
    def run_splitted_configs(run_function: Callable, configs: List[Dict]):
        app_config = AppConfig()
        with ProcessPoolExecutor(max_workers=app_config.THREADS) as executor:
            executor.map(run_function, configs)
