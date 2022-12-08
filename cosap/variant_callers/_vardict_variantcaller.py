import glob
import os
from pathlib import Path
from subprocess import PIPE, Popen, check_output, run
from typing import Dict, List, Union

from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class VarDictVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_vardict_command(
        cls, caller_config:Dict, library_paths:LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]
        tumor_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.TUMOR_SAMPLE_NAME
        ]

        command = [
            "vardict-java",
            "-G",
            library_paths.REF_FASTA,
            "-f 0.01",
            "-N",
            tumor_sample_name,
            "-b",
            f"'{tumor_bam}|{germline_bam}'",
            "|",
            "testsomatic.R",
        ]
        print(" ".join(command))
        return command

    @classmethod
    def _create_var2vcf_command(
        cls, caller_config:Dict, library_paths:LibraryPaths
    ) -> List:

        tumor_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.TUMOR_SAMPLE_NAME
        ]
        germline_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.GERMLINE_SAMPLE_NAME
        ]
        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]

        command = [
            "var2vcf_paired.pl",
            "-N",
            f"'{tumor_sample_name}|{germline_sample_name}'",
            "-f 0.05",
            ">",
            output_name,
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        vardict_command = cls._create_vardict_command(
            caller_config=caller_config, library_paths=library_paths
        )
        var2vcf_command = cls._create_var2vcf_command(
            caller_config=caller_config, library_paths=library_paths
        )

        vardict = Popen(vardict_command, stdout=PIPE, shell=True)
        var2vcf = check_output(var2vcf_command, stdin=vardict.stdout, shell=True)
        vardict.wait()
