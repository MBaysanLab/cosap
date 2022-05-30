from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import QualityControlKeys
from ._quality_controllers import _QualityControllable, _QualityController


class Qualimap(_QualityController, _QualityControllable):
    @classmethod
    def _create_qualimap_command(
        cls, qc_config=Dict, library_paths=LibraryPaths
    ) -> List:

        input_bam = qc_config[QualityControlKeys.INPUT]
        bed_file = (
            qc_config[QualityControlKeys.BED_FILE]
            if QualityControlKeys.BED_FILE in qc_config.keys()
            else None
        )
        raw_output = qc_config[QualityControlKeys.RAW_OUTPUT]
        output = qc_config[QualityControlKeys.OUTPUT]

        command = [
            "qualimap",
            "bamqc",
            "-bam",
            input_bam,
            "-outfile",
            output,
            "-outdir",
            raw_output,
        ]
        if bed_file is not None:
            command.extend(["--feature-file", bed_file])
        return command

    @classmethod
    def run_qualitycontroller(cls, qc_config=Dict, library_paths=LibraryPaths):
        qualimap_command = cls._create_qualimap_command(
            qc_config=qc_config, library_paths=library_paths
        )

        run(qualimap_command, cwd=qc_config[QualityControlKeys.OUTPUT_DIR])
