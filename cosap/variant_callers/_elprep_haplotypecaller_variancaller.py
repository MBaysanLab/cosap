import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class ElprepHaplotypecaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        command = [
            "elprep",
            "sfm",
            "--reference",
            library_paths.REF_ELFASTA,
            "--haplotypecaller",
            caller_config[VariantCallingKeys.GVCF_OUTPUT],
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        command = cls._create_run_command(
            caller_config=caller_config, library_paths=library_paths
        )

        run(command)
