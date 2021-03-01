import glob
import os
from abc import ABC, abstractmethod
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._sorting_config import SortingKeys
from .._utils import join_paths


class _Sorter(ABC):
    @abstractmethod
    def sort(self):
        pass


class _Sortable:
    pass


class SamtoolsSorter(_Sorter, _Sortable):
    @classmethod
    def _create_command(
        cls, sorting_config: Dict, app_config: AppConfig, library_paths: LibraryPaths
    ) -> str:
        command = [
            "samtools",
            "view",
            "-@",
            app_config.THREADS,
            "-bS",
            sorting_config[SortingKeys.INPUT],
            "|",
            "samtools",
            "sort",
            "-@",
            app_config.THREADS,
            "-o",
            sorting_config[SortingKeys.OUTPUT],
        ]
        return command

    @classmethod
    def sort(cls, sorting_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            sorting_config=sorting_config,
            app_config=app_config,
            library_paths=library_paths,
        )

        run(command, cwd=sorting_config.BAM_DIR)
