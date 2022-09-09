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
            "-Xmx16G",
            "HaplotypeCaller",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            germline_bam,
            "-O",
            output_name,
            "--native-pair-hmm-threads",
            str(AppConfig.THREADS),
            "-ERC",
            "GVCF"
        ]
        return command

    @classmethod
    def _create_cnnscorevariants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.FILTERED_VARIANTS_OUTPUT]
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        command = [
            "gatk",
            "CNNScoreVariants",
            "-I",
            input_bam,
            "-V",
            input_name,
            "-R",
            library_paths.REF_FASTA,
            "-O",
            output_name,
            "-tensor-type",
            "read-tensor",
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
            "--info-key",
            "CNN_1D",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_snp_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:

        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
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

        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
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

        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
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
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        haplotypecaller_command = cls._create_run_command(
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

        run(haplotypecaller_command)
        # run(get_snp_command)
        # run(get_indel_command)
        # run(get_other_variants_command)
