from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import ElprepKeys
from ._preprocessors import _PreProcessable, _Preprocessor


class ElprepPreprocess(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, elrep_config: Dict
    ) -> List:

        command = [
            "elprep",
            "filter",
            elrep_config[ElprepKeys.INPUT],
            elrep_config[ElprepKeys.OUTPUT],
            "--mark-duplicates",
            "--mark-optical-duplicates",
            f"{elrep_config[ElprepKeys.OUTPUT]}_metrics",
            "--bqsr",
            elrep_config[ElprepKeys.TABLE],
            "--reference",
            library_paths.REF_ELFASTA
        ]
        return command

    @classmethod
    def run_preprocessor(cls, elrep_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths, app_config=app_config, elrep_config=elrep_config
        )
        run(command)
