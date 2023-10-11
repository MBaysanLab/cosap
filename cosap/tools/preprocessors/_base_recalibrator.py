import os
from pathlib import Path
from subprocess import run
from typing import Dict

from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import BaseRecalibratorKeys
from ...pipeline_runner.runners import DockerRunner
from ._preprocessors import _PreProcessable, _Preprocessor


class BaseRecalibrator(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_table_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ):
        command = [
            "gatk",
            "BaseRecalibrator",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            calibration_config[BaseRecalibratorKeys.INPUT],
            "--known-sites",
            library_paths.MILLS_INDEL,
            "--known-sites",
            library_paths.DBSNP,
            "--known-sites",
            library_paths.ONE_THOUSAND_G,
            "-O",
            calibration_config[BaseRecalibratorKeys.TABLE],
        ]
        bed_file = (
            calibration_config[BaseRecalibratorKeys.BED_FILE]
            if BaseRecalibratorKeys.BED_FILE in calibration_config.keys()
            else None
        )
        if bed_file is not None:
            command.extend(["--intervals", bed_file])
        return command

    @classmethod
    def _create_calibration_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ):
        command = [
            "gatk",
            "ApplyBQSR",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            calibration_config[BaseRecalibratorKeys.INPUT],
            "--bqsr-recal-file",
            calibration_config[BaseRecalibratorKeys.TABLE],
            "-O",
            calibration_config[BaseRecalibratorKeys.OUTPUT],
        ]
        bed_file = (
            calibration_config[BaseRecalibratorKeys.BED_FILE]
            if BaseRecalibratorKeys.BED_FILE in calibration_config.keys()
            else None
        )
        if bed_file is not None:
            command.extend(["--intervals", bed_file])
        return command

    @classmethod
    def _create_table(cls, calibration_config: Dict, library_paths: LibraryPaths):
        command = cls._create_table_command(
            calibration_config=calibration_config, library_paths=library_paths
        )
        run(command)

    @classmethod
    def _apply_calibration(cls, calibration_config: Dict, library_paths: LibraryPaths):
        command = cls._create_calibration_command(
            calibration_config=calibration_config, library_paths=library_paths
        )
        run(command)

    @classmethod
    def _create_parabricks_bqsr_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ) -> list:

        command = [
            "pbrun",
            "bqsr",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            calibration_config[BaseRecalibratorKeys.INPUT],
            "--knownSites",
            library_paths.MILLS_INDEL,
            "--knownSites",
            library_paths.DBSNP,
            "--knownSites",
            library_paths.ONE_THOUSAND_G,
            "--out-recal-file",
            calibration_config[BaseRecalibratorKeys.TABLE],
        ]

        return command

    @classmethod
    def _create_parabricks_applybqsr_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ) -> list:

        command = [
            "pbrun",
            "applybqsr",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            calibration_config[BaseRecalibratorKeys.INPUT],
            "--in-recal-file",
            calibration_config[BaseRecalibratorKeys.TABLE],
            "--out-bam",
            calibration_config[BaseRecalibratorKeys.OUTPUT],
        ]

        return command

    @classmethod
    def run_preprocessor(cls, calibration_config: Dict, device="cpu"):
        library_paths = LibraryPaths()

        if device == "cpu":
            cls._create_table(
                calibration_config=calibration_config, library_paths=library_paths
            )
            cls._apply_calibration(
                calibration_config=calibration_config, library_paths=library_paths
            )
        elif device == "gpu":
            bqsr_command = cls._create_parabricks_bqsr_command(
                calibration_config=calibration_config, library_paths=library_paths
            )
            applybqsr_command = cls._create_parabricks_applybqsr_command(
                calibration_config=calibration_config, library_paths=library_paths
            )

            output_dir = os.path.abspath(
                os.path.dirname(calibration_config[BaseRecalibratorKeys.OUTPUT])
            )
            os.makedirs(output_dir, exist_ok=True)

            runner = DockerRunner()
            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(bqsr_command),
                workdir=str(Path(output_dir).parent),
            )
            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(applybqsr_command),
                workdir=str(Path(output_dir).parent),
            )
