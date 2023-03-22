import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MDUPKeys
from .._utils import join_paths
from ..memory_handler import MemoryHandler
from ._preprocessors import _PreProcessable, _Preprocessor


class MarkDuplicate(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(
        cls,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        mdup_config: Dict,
        memory_handler: MemoryHandler,
    ) -> List:

        germline_bam = memory_handler.get_path(mdup_config[MDUPKeys.INPUT])
        tmp_dir = memory_handler.get_temp_dir()

        command = [
            "gatk",
            "MarkDuplicatesSpark",
            "-I",
            germline_bam,
            "-O",
            mdup_config[MDUPKeys.OUTPUT],
            "-M",
            f"{mdup_config[MDUPKeys.OUTPUT]}_metrics",
            "--remove-all-duplicates",
            "--create-output-bam-index",
            "--spark-master",
            f"local[{app_config.MAX_THREADS_PER_JOB}]",
            "--tmp-dir",
            tmp_dir,
        ]
        return command

    @classmethod
    def run_preprocessor(cls, mdup_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        with MemoryHandler() as memory_handler:
            command = cls._create_command(
                library_paths=library_paths,
                app_config=app_config,
                mdup_config=mdup_config,
                memory_handler=memory_handler,
            )
            run(command)
