from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MDUPKeys


class MarkDuplicate:
    @classmethod
    def _create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, mdup_config: Dict
    ) -> List:

        command = [
            "java",
            f"-XX:ParallelGCThreads={app_config.THREADS}",
            "-jar",
            library_paths.PICARD,
            "MarkDuplicates",
            f"-I={mdup_config[MDUPKeys.INPUT]}",
            f"-O={mdup_config[MDUPKeys.OUTPUT]}",
            f"-M={mdup_config[MDUPKeys.METRICS]}",
            "REMOVE_DUPLICATES=true",
            "CREATE_INDEX=true",
        ]
        return command

    @classmethod
    def mark_duplicate(cls, mdup_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths, app_config=app_config, mdup_config=mdup_config
        )
        run(command, cwd=mdup_config[MDUPKeys.OUTPUT_DIR])