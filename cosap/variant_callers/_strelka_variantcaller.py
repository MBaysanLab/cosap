import gzip
import os
import shutil
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class StrelkaVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_strelka_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        command = [
            "configureStrelkaSomaticWorkflow.py",
            f"--tumorBam={tumor_bam}",
            f"--normalBam={germline_bam}",
            f"--referenceFasta={library_paths.REF_FASTA}",
        ]

        return command

    @classmethod
    def _create_run_strelka_workflow_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        command = [
            "StrelkaSomaticWorkflow/runWorkflow.py",
            "-m",
            "local",
            "-j",
            str(AppConfig.THREADS),
        ]
        return command

    @classmethod
    def _move_strelka_vcfs(cls, caller_config=Dict, library_paths=LibraryPaths) -> List:

        snvs = "StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"
        indels = "StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"

        snp_output_filename = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        indel_output_filename = caller_config[VariantCallingKeys.INDEL_OUTPUT]

        with gzip.open(snvs, "rb") as snv_in:
            with open(snp_output_filename, "wb") as snv_out:
                shutil.copyfileobj(snv_in, snv_out)

        with gzip.open(indels, "rb") as indel_in:
            with open(indel_output_filename, "wb") as indel_out:
                shutil.copyfileobj(indel_in, indel_out)

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

        cls._move_strelka_vcfs(caller_config=caller_config, library_paths=library_paths)
