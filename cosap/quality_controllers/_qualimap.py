from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import QualityControlKeys
from ._quality_controllers import _QualityControllable, _QualityController
import os

class Qualimap(_QualityController, _QualityControllable):
    @classmethod
    def _create_qualimap_command(
        cls, qc_config=Dict, app_config=AppConfig
    ) -> List:

        MAX_MEMORY_IN_GB = int(AppConfig.MAX_MEMORY_PER_JOBS // (1024.0**3))

        input_bam = qc_config[QualityControlKeys.INPUT]
        bed_file = (
            qc_config[QualityControlKeys.BED_FILE]
            if QualityControlKeys.BED_FILE in qc_config.keys()
            else None
        )
        raw_output = qc_config[QualityControlKeys.RAW_OUTPUT]
        output_file = qc_config[QualityControlKeys.OUTPUT]

        command = [
            "qualimap",
            "bamqc",
            "-bam",
            input_bam,
            f"--java-mem-size={MAX_MEMORY_IN_GB}G",
            "-outfile",
            output_file,
            "-outdir",
            raw_output,
            "-outformat",
            "PDF",
            "--output-genome-coverage",
            f"{raw_output}/coverage_histogram.txt",
            "--nt",
            str(app_config.MAX_THREADS_PER_JOB)

        ]
        if bed_file is not None:
            command.extend(["--feature-file", bed_file])
        return command

    @classmethod
    def run_qualitycontroller(cls, qc_config=Dict, library_paths=LibraryPaths):
        app_config = AppConfig()
        qualimap_command = cls._create_qualimap_command(
            qc_config=qc_config, app_config=app_config
        )

        run(qualimap_command)
