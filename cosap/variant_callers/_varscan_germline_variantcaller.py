import glob
import os
from pathlib import Path
from subprocess import PIPE, Popen, check_output, run
from typing import Dict, List, Union

from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class VarScanGermlineVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_samtools_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        command = [
            "samtools",
            "mpileup",
            "-B",
            "-q",
            "1",
            "-f",
            library_paths.REF_FASTA,
            germline_bam,
        ]
        return command

    @classmethod
    def _create_varscan_command(cls) -> List:

        command = [
            "varscan",
            "mpileup2snp",
            "--p-value",
            "99e-02"

        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        samtools_mpileup = cls._create_samtools_command(
            caller_config=caller_config, library_paths=library_paths
        )
        varscan_germline = cls._create_varscan_command()

        samtools = Popen(samtools_mpileup, stdout=PIPE)
        varscan = check_output(varscan_germline, stdin=samtools.stdout)
        samtools.wait()

        with open(caller_config[VariantCallingKeys.SNP_OUTPUT], "wb+") as vcf_file:
            vcf_file.write(varscan)
