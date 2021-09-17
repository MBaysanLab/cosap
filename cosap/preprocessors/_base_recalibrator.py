from typing import Dict
from subprocess import run
from .._library_paths import LibraryPaths
from .._pipeline_config import BaseRecalibratorKeys
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
            f"{library_paths.REF_DIR}/Homo_sapiens_assembly38.fasta",
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
        return command

    @classmethod
    def _create_calibration_command(
        cls, calibration_config: Dict, library_paths: LibraryPaths
    ):
        command = [
            "gatk",
            "ApplyBQSR",
            "-R",
            f"{library_paths.REF_DIR}/Homo_sapiens_assembly38.fasta",
            "-I",
            calibration_config[BaseRecalibratorKeys.INPUT],
            "--bqsr-recal-file",
            calibration_config[BaseRecalibratorKeys.TABLE],
            "-O",
            calibration_config[BaseRecalibratorKeys.OUTPUT],
        ]
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
    def run_preprocessor(cls, calibration_config: Dict):
        library_paths = LibraryPaths()

        cls._create_table(
            calibration_config=calibration_config, library_paths=library_paths
        )
        cls._apply_calibration(
            calibration_config=calibration_config, library_paths=library_paths
        )
