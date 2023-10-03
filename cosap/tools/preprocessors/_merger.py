from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MergingKeys


class BamMerger:
    @classmethod
    def _create_command(
        cls, merging_config: Dict, app_config: AppConfig, library_paths: LibraryPaths
    ) -> List:
        bam_files = [f"I={bam_file}" for bam_file in merging_config[MergingKeys.INPUTS]]
        command = [
            "java",
            f"-XX:ParallelGCThreads={app_config.MAX_THREADS_PER_JOB}",
            "-jar",
            library_paths.PICARD,
            "MergeSamFiles",
            *bam_files,
            f"O={merging_config[MergingKeys.OUTPUT]}",
        ]

        return command

    @classmethod
    def merge(cls, merging_config: Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            merging_config=merging_config,
            app_config=app_config,
            library_paths=library_paths,
        )

        run(command)
