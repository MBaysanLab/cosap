from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._sorting_config import MergingKeys

class BamMerger:
    @classmethod
    def _create_command(cls, merging_config:Dict, app_config: AppConfig, library_path:LibraryPaths) -> List:
        bam_files = merging_config[MergingKeys.Inputs]
        command = [
            "java",
            "-XX:ParallelGCThreads=",
            app_config.THREADS,
            "-jar",
            library_path.PICARD,
            "MergeSamFiles"]

        command+=[f"I={bam_file} " for bam_file in bam_files]
        command+=["O=",merging_config[MergingKeys.Output]]

        return command


    @classmethod
    def merge(cls, merging_config:Dict):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        command = cls._create_command(
            merging_config=merging_config,
            app_config=app_config,
            library_paths=library_paths
        )

        run(command, cwd=merging_config.BAM_DIR)