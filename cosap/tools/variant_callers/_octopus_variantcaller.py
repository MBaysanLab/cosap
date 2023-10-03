from subprocess import run
from typing import Dict

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class OctopusVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> list:
        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        germline_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.GERMLINE_SAMPLE_NAME
        ]

        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]

        command = [
            "octopus",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            germline_bam,
            tumor_bam,
            "--normal-sample",
            germline_sample_name,
            "-o",
            output_name,
            "--threads",
            str(app_config.MAX_THREADS_PER_JOB),
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        octopus_command = cls._create_run_command(
            caller_config=caller_config,
            library_paths=library_paths,
            app_config=app_config,
        )

        run(octopus_command)
