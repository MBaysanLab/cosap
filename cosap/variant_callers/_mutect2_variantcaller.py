import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class Mutect2VariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> list:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        germline_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.GERMLINE_SAMPLE_NAME
        ]
        tumor_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.TUMOR_SAMPLE_NAME
        ]

        output_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
        ]

        command = [
            "gatk",
            "Mutect2",
            "-R",
            library_paths.REF_DIR,
            "-I",
            germline_bam,
            "-tumor",
            tumor_sample_name,
            "-I",
            tumor_bam,
            "-normal",
            germline_sample_name,
            "-O",
            output_name,
        ]
        return command

    @classmethod
    def _create_get_snp_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> str:

        input_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
        ]
        output_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.SNP_OUTPUT
        ]

        command = [
            "gatk4",
            "SelectVariants",
            "-R",
            library_paths.REF_DIR,
            "-V",
            input_name,
            "--select-type-to-include",
            "SNP",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_indel_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> str:

        input_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
        ]
        output_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.INDEL_OUTPUT
        ]

        command = [
            "gatk",
            "SeleckVariants",
            "-R",
            library_paths.REF_DIR,
            "-V",
            input_name,
            "--select-type-to-include",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_other_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> str:

        input_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
        ]
        output_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.OTHER_VARIANTS_OUTPUT
        ]

        command = [
            "gatk",
            "SeleckVariants",
            "-R",
            library_paths.REF_DIR,
            "-V",
            input_name,
            "--select-type-to-exclude",
            "SNP",
            "--select-type-to-exclude",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        mutect_command = cls._create_run_command(
            caller_config=caller_config, library_paths=library_paths
        )
        get_snp_command = cls._create_get_snp_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )
        get_indel_command = cls._create_get_indel_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )
        get_other_variants_command = cls._create_get_other_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )

        run(mutect_command)
        run(get_snp_command)
        run(get_indel_command)
        run(get_other_variants_command)
