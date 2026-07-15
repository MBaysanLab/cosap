import os
from pathlib import Path
from subprocess import run
from typing import Dict

from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import BaseRecalibratorKeys
from ..._utils import convert_to_absolute_path
from ...runners.docker_runner import DockerRunner
from ._preprocessors import _PreProcessable, _Preprocessor


class BaseRecalibrator(_Preprocessor, _PreProcessable):
    @staticmethod
    def _create_output_path(calibration_config: Dict, output_key: str) -> str:
        output_path = calibration_config[output_key]
        if os.path.isabs(output_path):
            return output_path
        return os.path.join(
            calibration_config[BaseRecalibratorKeys.OUTPUT_DIR], output_path
        )

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
            convert_to_absolute_path(calibration_config[BaseRecalibratorKeys.INPUT]),
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
            convert_to_absolute_path(calibration_config[BaseRecalibratorKeys.BED_FILE])
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
            convert_to_absolute_path(calibration_config[BaseRecalibratorKeys.INPUT]),
            "--bqsr-recal-file",
            calibration_config[BaseRecalibratorKeys.TABLE],
            "-O",
            calibration_config[BaseRecalibratorKeys.OUTPUT],
            "--create-output-bam-index",
            "true",
        ]
        bed_file = (
            convert_to_absolute_path(calibration_config[BaseRecalibratorKeys.BED_FILE])
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
        run(
            command, cwd=calibration_config[BaseRecalibratorKeys.OUTPUT_DIR], check=True
        )

    @classmethod
    def _apply_calibration(cls, calibration_config: Dict, library_paths: LibraryPaths):
        command = cls._create_calibration_command(
            calibration_config=calibration_config, library_paths=library_paths
        )
        run(
            command, cwd=calibration_config[BaseRecalibratorKeys.OUTPUT_DIR], check=True
        )

    @classmethod
    def _create_parabricks_bqsr_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ) -> list:
        recal_table = cls._create_output_path(
            calibration_config, BaseRecalibratorKeys.TABLE
        )
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
            recal_table,
        ]

        return command

    @classmethod
    def _create_parabricks_applybqsr_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ) -> list:
        recal_table = cls._create_output_path(
            calibration_config, BaseRecalibratorKeys.TABLE
        )
        output_bam = cls._create_output_path(
            calibration_config, BaseRecalibratorKeys.OUTPUT
        )
        command = [
            "pbrun",
            "applybqsr",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            calibration_config[BaseRecalibratorKeys.INPUT],
            "--in-recal-file",
            recal_table,
            "--out-bam",
            output_bam,
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

            output_dir = calibration_config[BaseRecalibratorKeys.OUTPUT_DIR]
            if not os.path.isabs(output_dir):
                workdir = str(
                    Path(calibration_config[BaseRecalibratorKeys.WORKDIR]).absolute()
                )
                output_dir = os.path.join(
                    workdir,
                    output_dir,
                )
            else:
                workdir = str(
                    Path(calibration_config[BaseRecalibratorKeys.WORKDIR]).absolute()
                )
            os.makedirs(output_dir, exist_ok=True)

            runner = DockerRunner(device=device)
            return_code = runner.run(
                image=DockerImages.PARABRICKS,
                command=bqsr_command,
                workdir=workdir,
            )
            if return_code not in (0, None):
                raise RuntimeError(
                    f"Parabricks BQSR failed with return code {return_code}"
                )
            return_code = runner.run(
                image=DockerImages.PARABRICKS,
                command=applybqsr_command,
                workdir=workdir,
            )
            if return_code not in (0, None):
                raise RuntimeError(
                    f"Parabricks ApplyBQSR failed with return code {return_code}"
                )
