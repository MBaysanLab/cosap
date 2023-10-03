from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class MuseVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_muse_output_prefix(cls, output_name: str):
        return output_name.strip(".vcf")

    @classmethod
    def _create_muse_call_command(
        cls, caller_config: Dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> List:
        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]
        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_prefix = cls._create_muse_output_prefix(output_name)

        command = [
            "MuSE",
            "call",
            "-O",
            output_prefix,
            "-f",
            library_paths.REF_FASTA,
            tumor_bam,
            germline_bam,
            "-n",
            str(app_config.MAX_THREADS_PER_JOB),
        ]
        return command

    @classmethod
    def _create_muse_sump_command(
        cls, caller_config: Dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> List:
        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_prefix = cls._create_muse_output_prefix(output_name)

        command = [
            "MuSE",
            "sump",
            "-I",
            f"{output_prefix}.MuSE.txt",
            "-E",
            "-O",
            output_name,
            "-D",
            f"{library_paths.DBSNP}.gz",
            "-n",
            str(app_config.MAX_THREADS_PER_JOB),
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        call_command = cls._create_muse_call_command(
            caller_config=caller_config,
            library_paths=library_paths,
            app_config=app_config,
        )
        sump_command = cls._create_muse_sump_command(
            caller_config=caller_config,
            library_paths=library_paths,
            app_config=app_config,
        )

        run(call_command)
        run(sump_command)
