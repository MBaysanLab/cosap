import gzip
import os
import shutil
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller


class MantaVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_manta_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        command = [
            "configManta.py",
            f"--tumorBam={tumor_bam}",
            f"--referenceFasta={library_paths.REF_FASTA}",
        ]
        bed_file = (
            caller_config[VariantCallingKeys.BED_FILE]
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )
        if bed_file is not None:
            command.extend(["--callRegions", bed_file, "--exome"])

        return command

    @classmethod
    def _create_run_manta_workflow_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:

        command = [
            "MantaWorkflow/runWorkflow.py",
            "-m",
            "local",
            "-j",
            str(AppConfig.MAX_THREADS_PER_JOB),
        ]
        return command

    @classmethod
    def _move_manta_vcfs(cls, caller_config=Dict, library_paths=LibraryPaths) -> List:

        tumor_sv = "MantaWorkflow/results/variants/tumorSV.vcf.gz"

        sv_output_filename = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]

        with gzip.open(tumor_sv, "rb") as snv_in:
            with open(sv_output_filename, "wb") as snv_out:
                shutil.copyfileobj(snv_in, snv_out)

    @classmethod
    def call_variants(cls, caller_config=Dict):
        library_paths = LibraryPaths()

        manta_command = cls._create_manta_command(
            caller_config=caller_config, library_paths=library_paths
        )
        manta_run_wf_command = cls._create_run_manta_workflow_command(
            caller_config=caller_config, library_paths=library_paths
        )

        run(manta_command)
        run(manta_run_wf_command)

        cls._move_manta_vcfs(caller_config=caller_config, library_paths=library_paths)
