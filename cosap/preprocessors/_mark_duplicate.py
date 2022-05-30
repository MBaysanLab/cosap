from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MDUPKeys
from ._preprocessors import _PreProcessable, _Preprocessor


class MarkDuplicate(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, mdup_config: Dict
    ) -> List:

        command = [
            "picard",
            "MarkDuplicates",
            "--INPUT",
            mdup_config[MDUPKeys.INPUT],
            "--OUTPUT",
            mdup_config[MDUPKeys.OUTPUT],
            "--METRICS_FILE",
            f"{mdup_config[MDUPKeys.OUTPUT]}_metrics",
            "--REMOVE_DUPLICATES",
            "true",
            "--CREATE_INDEX",
            "true",
        ]
        return command

    @classmethod
    def run_preprocessor(cls, mdup_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths, app_config=app_config, mdup_config=mdup_config
        )
        run(command, cwd=mdup_config[MDUPKeys.OUTPUT_DIR])
