import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class HaplotypeCallerVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        output_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "--java-options",
            "'-Xmx16G'",
            "HaplotypeCaller",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            germline_bam,
            "-O",
            output_name,
            "--native-pair-hmm-threads",
            str(AppConfig.THREADS),
        ]
        return command

    @classmethod
    def _create_cnnscorevariants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.FILTERED_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "CNNScoreVariants",
            "-V",
            input_name,
            "-R",
            library_paths.REF_FASTA,
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_filter_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        input_name = caller_config[VariantCallingKeys.FILTERED_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.FILTERED_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "FilterVariantTranches",
            "-V",
            input_name,
            "--resource",
            library_paths.MILLS_INDEL,
            "--resource",
            library_paths.DBSNP,
            "--resource",
            library_paths.ONE_THOUSAND_G,
            "--info-key CNN_1D",
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
        cnn_score_command = cls._create_cnnscorevariants_command(
            caller_config=caller_config, library_paths=library_paths
        )
        filter_command = cls._create_filter_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )

        run(mutect_command)
        run(cnn_score_command)
        run(filter_command)
