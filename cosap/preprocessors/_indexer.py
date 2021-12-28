from pathlib import Path
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import IndexingKeys


class BamIndexer:
    @classmethod
    def _create_command(cls, library_paths: LibraryPaths, indexing_config: Dict) -> List:
        # TODO: give output from config
        command = [
            "java",
            "-jar",
            library_paths.PICARD,
            "BuildBamIndex",
<<<<<<< HEAD
            "-I=",
            indexing_config[IndexingKeys.INPUT],
            "-O",
            indexing_config[IndexingKeys.OUTPUT]
=======
            f"-I={bam_file}",
>>>>>>> e3fefa353b8983d2a7daaa176e98dc4e817b4aae
        ]
        return command

    @classmethod
    def _check_is_sorted(cls, indexing_config: Dict) -> List:
        command = ["samtools", "stats", indexing_config[IndexingKeys.INPUT], "|", "grep", "is sorted:"]
        return command

    @classmethod
    def create_index(cls, indexing_config: Dict):
        library_paths = LibraryPaths()

<<<<<<< HEAD
=======
        bam_file = indexing_config[IndexingKeys.INPUT]
>>>>>>> e3fefa353b8983d2a7daaa176e98dc4e817b4aae

        check_is_sorted_command = cls._check_is_sorted(indexing_config=indexing_config)
        # TODO: exception handling
        is_sorted = run(check_is_sorted_command).stdout
        if not is_sorted == "1":
            raise Exception("BAM file must be sorted before indexing")

<<<<<<< HEAD
        command = cls._create_command(library_paths=library_paths, indexing_config=indexing_config)
        run(command, cwd=indexing_config.BAM_DIR)
=======
        command = cls._create_command(library_paths=library_paths, bam_file=bam_file)
        run(command, cwd=indexing_config[IndexingKeys.BAM_DIR])
>>>>>>> e3fefa353b8983d2a7daaa176e98dc4e817b4aae
