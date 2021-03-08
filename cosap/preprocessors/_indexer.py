import glob
import os
from pathlib import Path
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import SortingKeys


class BamIndexer:
    @classmethod
    def _create_command(cls, library_paths: LibraryPaths, indexing_config: Dict) -> List:
        # TODO: give output from config
        command = [
            "java",
            "-jar",
            library_paths.PICARD,
            "BuildBamIndex",
            "-I=",
            indexing_config[IndexingKeys.INPUT],
            "-O",
            indexing_config[IndexingKeys.OUTPUT]
        ]
        return command

    @classmethod
    def _check_is_sorted(cls, indexing_config: Dict) -> List:
        command = ["samtools", "stats", indexing_config[IndexingKeys.INPUT], "|", "grep", "is sorted:"]
        return command

    # TODO: Preprocessing config instead of sorting config?
    @classmethod
    def create_index(cls, indexing_config: Dict):
        library_paths = LibraryPaths()


        check_is_sorted_command = cls._check_is_sorted(indexing_config=indexing_config)
        # TODO: exception handling
        is_sorted = run(check_is_sorted_command).stdout
        if not is_sorted == "1":
            raise Exception("BAM file must be sorted before indexing")

        command = cls._create_command(library_paths=library_paths, indexing_config=indexing_config)
        run(command, cwd=indexing_config.BAM_DIR)
