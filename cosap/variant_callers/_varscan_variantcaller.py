import glob
import os
from pathlib import Path
from subprocess import PIPE, Popen, check_output, run
from typing import Dict, List, Union

from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class VarScanVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_samtools_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        command = [
            "samtools",
            "mpileup",
            "-B",
            "-q",
            "1",
            "-f",
            library_paths.REF_FASTA,
            germline_bam,
            tumor_bam,
        ]
        return command

    @classmethod
    def _create_varscan_somatic_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        snp_output_name = caller_config[VariantCallingKeys.SNP_OUTPUT]
        indel_output_name = caller_config[VariantCallingKeys.INDEL_OUTPUT]

        command = [
            "varscan",
            "somatic",
            "--output-snp",
            snp_output_name,
            "--output-indel",
            indel_output_name,
            "--mpileup",
            "1",
            "--min-coverage-normal 8",
            "--min-coverage-tumor 6",
            "--min-var-freq 0.10",
            "--min-freq-for-hom 0.75",
            "--normal-purity 1.0",
            "--tumor-purity 1.00",
            "--p-value 0.99",
            "--somatic-p-value 0.05",
            "--strand-filter 0",
            "--output-vcf",
        ]
        return command

    @classmethod
    def _create_process_somatic_command(
        cls,
        vcf=Union[str, Path],
    ) -> List:
        command = [
            "varscan",
            "processSomatic",
            vcf,
            "--min-tumor-freq 0.10",
            "--max-normal-freq 0.05",
            "--p-value 0.07",
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        samtools_mpileup = cls._create_samtools_command(
            caller_config=caller_config, library_paths=library_paths
        )
        varscan_somatic = cls._create_varscan_somatic_command(
            caller_config=caller_config, library_paths=library_paths
        )

        samtools = Popen(samtools_mpileup, stdout=PIPE)
        varscan = check_output(varscan_somatic, stdin=samtools.stdout)
        samtools.wait()

        unfiltered_vcfs = [
            caller_config[VariantCallingKeys.SNP_OUTPUT],
            caller_config[VariantCallingKeys.INDEL_OUTPUT],
        ]

        for vcf_file in unfiltered_vcfs:
            process_somatic_command = cls._create_process_somatic_command(
                vcf=vcf_file,
            )
            run(process_somatic_command)
