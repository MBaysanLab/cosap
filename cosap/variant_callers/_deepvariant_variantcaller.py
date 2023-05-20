import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class DeepVariantVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def create_run_deepvariant_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        output_vcf = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_gvcf = caller_config[VariantCallingKeys.GVCF_OUTPUT]

        command = [
            "run_deepvariant",
            "--model_type=WGS",
            f"--ref={library_paths.REF_FASTA}",
            f"--reads={input_bam}",
            f"--output_vcf={output_vcf}",
            f"--output_gvcf={output_gvcf}",
            f"--num_shards={AppConfig.MAX_THREADS_PER_JOB}",
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()
        command = cls.create_run_deepvariant_command(caller_config, library_paths)
        run(command)