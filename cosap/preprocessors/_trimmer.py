from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import TrimmingKeys
from ._preprocessors import _Preprocessor, _PreProcessable


class Trimmer(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, trimmer_config: Dict
    ) -> List:

        fastq_inputs = " ".join(
            [fastq for fastq in trimmer_config[TrimmingKeys.INPUT].values()]
        )
        command = [
            "fastp",
            "-w",
            str(app_config.THREADS),
            "--in1",
            fastq_inputs[0],
            "--in2",
            fastq_inputs[1],
            "--out1",
            trimmer_config[TrimmingKeys.OUTPUT]["1"],
            "--out2",
            trimmer_config[TrimmingKeys.OUTPUT]["2"],
        ]
        return command

    @classmethod
    def run_preprocessor(cls, trimmer_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths, app_config=app_config, trimmer_config=trimmer_config
        )
        run(command)
