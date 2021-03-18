from pathlib import Path
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import IndexingKeys


class BamIndexer:
    @classmethod
    def _create_command(cls, library_paths: LibraryPaths, bam_file: Path) -> List:
        # TODO: give output from config
        command = [
            "java",
            "-jar",
            library_paths.PICARD,
            "BuildBamIndex",
            f"-I={bam_file}",
        ]
        return command

    @classmethod
    def _check_is_sorted(cls, bam_file: Path) -> List:
        command = ["samtools", "stats", bam_file, "|", "grep", "is sorted:"]
        return command

    @classmethod
    def create_index(cls, indexing_config: Dict):
        library_paths = LibraryPaths()

        bam_file = indexing_config[IndexingKeys.INPUT]

        check_is_sorted_command = cls._check_is_sorted(bam_file=bam_file)
        # TODO: exception handling
        is_sorted = run(check_is_sorted_command).stdout
        if not is_sorted == "1":
            raise Exception("BAM file must be sorted before indexing")

        command = cls._create_command(library_paths=library_paths, bam_file=bam_file)
        run(command, cwd=indexing_config[IndexingKeys.BAM_DIR])
