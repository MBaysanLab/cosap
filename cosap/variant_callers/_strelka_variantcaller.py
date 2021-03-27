import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class Strelka2VariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_strelka_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:
        bam_paths = cls._get_bam_paths(caller_config)
        command = [
            library_paths.STRELKA,
            "--normalBam",
            bam_paths["germline_bam_path"],
            "--tumorBam",
            bam_paths["tumor_bam_path"],
            "--referenceFasta",
            library_paths.REF_DIR,
            "--runDir",
            caller_config.VCF_OUTPUT_DIR,
            "--exome",
            "--disableEVS",
        ]

        return command

    @classmethod
    def call_variants(cls, caller_config=Dict):
        library_paths = LibraryPaths()
        strelka_command = cls._create_strelka_command(
            caller_config=caller_config, library_paths=library_paths
        )

        run(strelka_command, cwd=caller_config[VariantCallingKeys.OUTPUT_DIR])
