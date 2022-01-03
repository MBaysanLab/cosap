import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller
from .._config import AppConfig


class StrelkaVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_strelka_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        command = [
            "configureStrelkaSomaticWorkflow.py",
            "--tumorBam=",
            tumor_bam,
            "--normalBam=",
            germline_bam,
            "--referenceFasta=",
            library_paths.REF_DIR,
        ]

        return command

    @classmethod
    def _create_run_strelka_workflow_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        command = [
            "runWorkflow.py",
            "-m",
            "local",
            "-j",
            AppConfig.THREADS
        ]
        return command


    @classmethod
    def call_variants(cls, caller_config=Dict):
        library_paths = LibraryPaths()

        strelka_command = cls._create_strelka_command(
            caller_config=caller_config, library_paths=library_paths
        )
        strelka_run_wf_command = cls._create_run_strelka_workflow_command(
            caller_config=caller_config, library_paths=library_paths
        )
        run(strelka_command)
        run(strelka_run_wf_command)
