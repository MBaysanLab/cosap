import gzip
import os
import shutil
from pathlib import Path

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ..._utils import join_paths
from ...memory_handler import MemoryHandler
from ...runners.runners import DockerRunner
from ._variantcallers import _Callable, _VariantCaller


class MantaVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_manta_command(
        cls,
        caller_config: dict,
        library_paths: LibraryPaths,
    ) -> tuple[str, list]:
        tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]
        run_dir = caller_config[VariantCallingKeys.OUTPUT_DIR]

        command = [
            "configManta.py",
            f"--tumorBam={tumor_bam}",
            f"--referenceFasta={library_paths.REF_FASTA}",
            f"--runDir={run_dir}",
        ]
        bed_file = (
            caller_config[VariantCallingKeys.BED_FILE]
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )
        if bed_file is not None:
            command.extend(["--callRegions", bed_file, "--exome"])

        return run_dir, command

    @classmethod
    def _create_run_manta_workflow_command(
        cls, caller_config: dict, library_paths: LibraryPaths, rundir=str
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
    def _move_manta_vcfs(
        cls, caller_config: dict, library_paths: LibraryPaths, rundir: str
    ) -> list:
        tumor_sv = f"{rundir}/results/variants/tumorSV.vcf.gz"

        sv_output_filename = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT],
        )

        with gzip.open(tumor_sv, "rb") as snv_in:
            with open(sv_output_filename, "wb") as snv_out:
                shutil.copyfileobj(snv_in, snv_out)

    @classmethod
    def call_variants(cls, caller_config: dict, device: str = "cpu"):
        library_paths = LibraryPaths()

        workdir = os.getcwd()

        rundir, manta_command = cls._create_manta_command(
            caller_config=caller_config, library_paths=library_paths
        )
        manta_run_wf_command = cls._create_run_manta_workflow_command(
            caller_config=caller_config, library_paths=library_paths, rundir=rundir
        )

        docker_runner = DockerRunner(device=device)
        docker_runner.run(
            DockerImages.MANTA,
            " ".join(manta_command),
            workdir=workdir,
        )
        docker_runner.run(
            DockerImages.MANTA,
            " ".join(manta_run_wf_command),
            workdir=workdir,
        )
        cls._move_manta_vcfs(
            caller_config=caller_config,
            library_paths=library_paths,
            rundir=rundir,
        )
