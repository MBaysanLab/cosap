from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import QualityControlKeys
from ._quality_controllers import _QualityControllable, _QualityController


class Mosdepth(_QualityController, _QualityControllable):
    @classmethod
    def _create_mosdepth_command(
        cls, qc_config=Dict, library_paths=LibraryPaths
    ) -> List:

        input_bam = qc_config[QualityControlKeys.INPUT]
        bed_file = (
            qc_config[QualityControlKeys.BED_FILE]
            if QualityControlKeys.BED_FILE in qc_config.keys()
            else None
        )
        output = qc_config[QualityControlKeys.OUTPUT]

        command = [
            "mosdepth",
            "-n",
            "--fast-mode",
            "-t",
            "4",
            output,
            input_bam,
        ]
        if bed_file:
            command.extend(
                [
                    "--by",
                    bed_file
                ]
            )

        return command

    @classmethod
    def run_qualitycontroller(cls, qc_config=Dict, library_paths=LibraryPaths):
        mosdepth_command = cls._create_mosdepth_command(
            qc_config=qc_config, library_paths=library_paths
        )

        run(mosdepth_command)
