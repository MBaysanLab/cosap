from subprocess import run
from typing import Dict, List

from ..._pipeline_config import IndexingKeys


class BamIndexer:
    @classmethod
    def _create_command(cls, indexing_config: Dict) -> List:
        # TODO: give output from config
        command = [
            "picard",
            "BuildBamIndex",
            "--INPUT",
            indexing_config[IndexingKeys.INPUT],
            "--OUTPUT",
            indexing_config[IndexingKeys.OUTPUT],
        ]
        return command

    @classmethod
    def run_preprocessor(cls, indexing_config: Dict):
        command = cls._create_command(indexing_config=indexing_config)
        run(command)
