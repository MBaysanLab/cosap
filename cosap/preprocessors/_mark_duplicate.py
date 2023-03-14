from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MDUPKeys
from ._preprocessors import _PreProcessable, _Preprocessor
from .._utils import join_paths
import os

class MarkDuplicate(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, mdup_config: Dict
    ) -> List:
        
        tmpdir_dir = join_paths(
            os.path.dirname(mdup_config[MDUPKeys.OUTPUT]),
            "tmp")

        os.makedirs(tmpdir_dir, exist_ok=True)

        command = [
            "gatk",
            "MarkDuplicatesSpark",
            "-I",
            mdup_config[MDUPKeys.INPUT],
            "-O",
            mdup_config[MDUPKeys.OUTPUT],
            "-M",
            f"{mdup_config[MDUPKeys.OUTPUT]}_metrics",
            "--remove-all-duplicates",
            "--create-output-bam-index",
            "--spark-master",
            f"local[{app_config.MAX_THREADS_PER_JOB}]",
            "--tmp-dir",
            tmpdir_dir
        ]
        return command

    @classmethod
    def run_preprocessor(cls, mdup_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths, app_config=app_config, mdup_config=mdup_config
        )
        run(command)
