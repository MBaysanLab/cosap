from subprocess import run
from typing import Dict, List

from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class SomaticSniperVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_somaticSniper_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:
        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        command = [
            "bam-somaticsniper",
            "-G",
            "-L",
            "-f",
            library_paths.REF_FASTA,
            "-F",
            "vcf",
            tumor_bam,
            germline_bam,
            output_name,
        ]

        return command

    @classmethod
    def call_variants(cls, caller_config=Dict, device: str = "cpu"):
        library_paths = LibraryPaths()

        somatic_sniper_command = cls._create_somaticSniper_command(
            caller_config=caller_config, library_paths=library_paths
        )
        run(somatic_sniper_command, cwd=caller_config[VariantCallingKeys.OUTPUT_DIR])
