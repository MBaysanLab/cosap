import os
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import TrimmingKeys
from ._preprocessors import _PreProcessable, _Preprocessor


class Trimmer(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, trimmer_config: Dict
    ) -> List:
        fastq_inputs = [fastq for fastq in trimmer_config[TrimmingKeys.INPUT].values()]
        command = [
            "fastp",
            "-w",
            str(app_config.MAX_THREADS_PER_JOB),
            "--in1",
            fastq_inputs[0],
            "--in2",
            fastq_inputs[1],
            "--out1",
            trimmer_config[TrimmingKeys.OUTPUT]["1"],
            "--out2",
            trimmer_config[TrimmingKeys.OUTPUT]["2"],
            "-h",
            os.devnull,
            "-j",
            trimmer_config[TrimmingKeys.REPORT_OUTPUT],
        ]
        return command

    @classmethod
    def run_preprocessor(cls, trimmer_config: Dict, *args, **kwargs):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths,
            app_config=app_config,
            trimmer_config=trimmer_config,
        )
        run(command)
