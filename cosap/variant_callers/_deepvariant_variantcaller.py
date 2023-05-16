import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class DeepVariantVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> list:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        output_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]

        command = [
            "python",
            library_paths.DEEPVARIANT_RUNNER,
            "--model_type=WGS",
            f"--ref={library_paths.REF_FASTA}",
            f"--reads={germline_bam}",
            f"--output_vcf={output_name}",
            f"--output_gvcf={output_name}.gvcf",
            f"--num_shards={app_config.MAX_THREADS_PER_JOB}",
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        deepvariant_command = cls._create_run_command(
            caller_config=caller_config,
            library_paths=library_paths,
            app_config=app_config,
        )

        run(deepvariant_command)
