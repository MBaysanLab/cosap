import glob
import os
from abc import ABC, abstractmethod
from subprocess import run
from typing import Dict, List
from pathlib import Path

from .._library_paths import LibraryPaths


class _Indexer(ABC):
    @abstractmethod
    def create_index():
        pass


class _Indexable(ABC):
    pass


class BamIndexer(_Indexer, _Indexable):
    @classmethod
    def _create_command(cls, library_paths: LibraryPaths, bam_file: Path) -> List:
        command = [
            "java",
            "-jar",
            library_paths.PICARD,
            "BuildBamIndex",
            "-I=",
            bam_file,
        ]
        return command

    @classmethod
    def _check_is_sorted(cls, bam_file: Path) -> List:
        command = ["samtools", "stats", bam_file, "|", "grep", "is sorted:"]

        return command

    # TODO: Preprocessing config instead of sorting config?
    @classmethod
    def create_index(cls, sorting_config: Dict, bam_file: Path):
        library_paths = LibraryPaths()

        check_is_sorted_command = cls._check_is_sorted(bam_file=bam_file)
        is_sorted = run(check_is_sorted_command).stdout
        if not is_sorted == "1":
            raise Exception("BAM file must be sorted before indexing")

        command = cls._create_command(library_paths=library_paths, bam_file=bam_file)
        run(command, cwd=sorting_config.BAM_DIR)
