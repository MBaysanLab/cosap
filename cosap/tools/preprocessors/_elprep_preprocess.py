from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import ElprepKeys
from ._preprocessors import _PreProcessable, _Preprocessor


class ElprepPreprocess(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(cls, library_paths: LibraryPaths, elprep_config: Dict) -> List:
        command = [
            "elprep",
            "sfm",
            elprep_config[ElprepKeys.INPUT],
            elprep_config[ElprepKeys.OUTPUT],
            "--mark-duplicates",
            "--mark-optical-duplicates",
            f"{elprep_config[ElprepKeys.OUTPUT]}_metrics",
            "--bqsr",
            elprep_config[ElprepKeys.TABLE],
            "--reference",
            library_paths.REF_ELFASTA,
        ]
        return command

    @classmethod
    def _index_output(cls, elprep_config: Dict) -> List:
        command = [
            "samtools",
            "index",
            elprep_config[ElprepKeys.OUTPUT],
        ]
        return command

    @classmethod
    def run_preprocessor(cls, elprep_config: Dict):
        library_paths = LibraryPaths()

        command = cls._create_command(
            library_paths=library_paths,
            elprep_config=elprep_config,
        )
        index_command = cls._index_output(elprep_config=elprep_config)
        run(command)
        run(index_command)
