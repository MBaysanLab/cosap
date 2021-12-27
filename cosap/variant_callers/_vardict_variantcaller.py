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
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]
        tumor_sample_name = caller_config[VariantCallingKeys.PARAMS][
                VariantCallingKeys.TUMOR_SAMPLE_NAME
            ]

        output_name = caller_config[VariantCallingKeys.SNP_OUTPUT]
        command = [
            "VarDict",
            "-G",
            LibraryPaths.REF_FASTA,
            "-N",
            tumor_sample_name,
            "-b",
            f"{tumor_bam}|{germline_bam}"
        ]