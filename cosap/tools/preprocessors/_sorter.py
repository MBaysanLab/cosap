from subprocess import run
from typing import Dict

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import SortingKeys
from ..._utils import join_paths


class SamtoolsSorter:
    @classmethod
    def _create_command(
        cls, sorting_config: Dict, app_config: AppConfig, library_paths: LibraryPaths
    ) -> str:
        command = [
            "samtools",
            "sort",
            sorting_config[SortingKeys.INPUT],
            "-@",
            str(app_config.MAX_THREADS_PER_JOB),
            "-o",
            sorting_config[SortingKeys.OUTPUT],
        ]
        if sorting_config[SortingKeys.SORTING_METHOD] == "name":
            command.append("-n")
        return command

    @classmethod
    def run_preprocessor(cls, sorting_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            sorting_config=sorting_config,
            app_config=app_config,
            library_paths=library_paths,
        )

        run(command)
