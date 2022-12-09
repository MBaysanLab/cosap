import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class Mutect2VariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        MAX_MEMORY_IN_GB = int(AppConfig.MAX_MEMORY_PER_JOBS // (1024**3))

        germline_bam = (
            caller_config[VariantCallingKeys.GERMLINE_INPUT]
            if VariantCallingKeys.GERMLINE_INPUT in caller_config.keys()
            else None
        )
        tumor_bam = (
            caller_config[VariantCallingKeys.TUMOR_INPUT]
            if VariantCallingKeys.TUMOR_INPUT in caller_config.keys()
            else None
        )

        output_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]



        command = [
            "gatk",
            "--java-options",
            f"-Xmx{MAX_MEMORY_IN_GB}G",
            "Mutect2",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            tumor_bam,
            "-O",
            output_name,
        ]

        if germline_bam is not None:
            germline_sample_name = caller_config[VariantCallingKeys.PARAMS][
                VariantCallingKeys.GERMLINE_SAMPLE_NAME
            ]
            command.extend(
                [
                    "-I",
                    germline_bam,
                    "-normal",
                    germline_sample_name,
                ]
            )
        return command

    @classmethod
    def _create_get_snp_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        input_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.SNP_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
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
    ) -> List:

        input_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.INDEL_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
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
    ) -> List:

        input_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.OTHER_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
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
    def _filter_mutect_calls(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "FilterMutectCalls",
            "-V",
            input_name,
            "-R",
            library_paths.REF_FASTA,
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
        filter_command = cls._filter_mutect_calls(
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
        run(filter_command)
        run(get_snp_command)
        run(get_indel_command)
        run(get_other_variants_command)
