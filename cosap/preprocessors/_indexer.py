from pathlib import Path
from subprocess import run
from typing import Dict, List

from .._pipeline_config import IndexingKeys


class BamIndexer:
    @classmethod
    def _create_command(cls, indexing_config: Dict) -> List:
        # TODO: give output from config
        command = [
            "picard" "BuildBamIndex",
            "-I=",
            indexing_config[IndexingKeys.INPUT],
            "-O",
            indexing_config[IndexingKeys.OUTPUT],
        ]
        return command

    @classmethod
    def _check_is_sorted(cls, indexing_config: Dict) -> List:
        command = [
            "samtools",
            "stats",
            indexing_config[IndexingKeys.INPUT],
            "|",
            "grep",
            "is sorted:",
        ]
        return command

    @classmethod
    def create_index(cls, indexing_config: Dict):

        check_is_sorted_command = cls._check_is_sorted(indexing_config=indexing_config)
        # TODO: exception handling
        is_sorted = run(check_is_sorted_command).stdout
        if not is_sorted == "1":
            raise Exception("BAM file must be sorted before indexing")

        command = cls._create_command(indexing_config=indexing_config)
        run(command, cwd=indexing_config.BAM_DIR)
