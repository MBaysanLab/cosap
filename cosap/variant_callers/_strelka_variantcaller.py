import gzip
import os
import shutil
from subprocess import run

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller
from ..memory_handler import MemoryHandler


class StrelkaVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_strelka_command(
        cls, caller_config:dict, library_paths:LibraryPaths, memory_handler:MemoryHandler
    ) -> tuple[str,list]:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]
        tmp_dir = memory_handler.get_temp_dir()

        command = [
            "configureStrelkaSomaticWorkflow.py",
            f"--tumorBam={tumor_bam}",
            f"--normalBam={germline_bam}",
            f"--referenceFasta={library_paths.REF_FASTA}",
            f"--runDir={tmp_dir}",
        ]

        return tmp_dir, command

    @classmethod
    def _create_run_strelka_workflow_command(
        cls, caller_config:dict, library_paths:LibraryPaths, rundir=str
    ) -> list:

        command = [
            f"{rundir}/runWorkflow.py",
            "-m",
            "local",
            "-j",
            str(AppConfig.MAX_THREADS_PER_JOB),
        ]
        return command

    @classmethod
    def _move_strelka_vcfs(cls, caller_config:dict, library_paths:LibraryPaths, rundir:str) -> list:

        snvs = f"{rundir}/results/variants/somatic.snvs.vcf.gz"
        indels = f"{rundir}/results/variants/somatic.indels.vcf.gz"

        snp_output_filename = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        indel_output_filename = caller_config[VariantCallingKeys.INDEL_OUTPUT]

        with gzip.open(snvs, "rb") as snv_in:
            with open(snp_output_filename, "wb") as snv_out:
                shutil.copyfileobj(snv_in, snv_out)

        with gzip.open(indels, "rb") as indel_in:
            with open(indel_output_filename, "wb") as indel_out:
                shutil.copyfileobj(indel_in, indel_out)

    @classmethod
    def call_variants(cls, caller_config:dict):
        library_paths = LibraryPaths()

        with MemoryHandler() as memory_handler:
            rundir, strelka_command = cls._create_strelka_command(
                caller_config=caller_config, library_paths=library_paths
            )
            strelka_run_wf_command = cls._create_run_strelka_workflow_command(
                caller_config=caller_config, library_paths=library_paths, rundir=rundir
            )

            run(strelka_command)
            run(strelka_run_wf_command)

            cls._move_strelka_vcfs(caller_config=caller_config, library_paths=library_paths, rundir=rundir)
