from pathlib import Path
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import SplitKeys


class SplitbyCHR:
    @classmethod
    def _create_command(
        cls,
        library_paths: LibraryPaths,
        bam_filename: Path,
        chromosome: str,
        output_filename: str,
    ) -> List:
        command = [
            "samtools",
            "view",
            "-bh",
            bam_filename,
            f"chr{chromosome}",
            ">",
            output_filename,
        ]
        return command

    @classmethod
    def split_by_chr(cls, split_config: Dict):
        library_paths = LibraryPaths()
        bam_files = split_config[SplitKeys.INPUT]
        output_files = split_config[SplitKeys.OUTPUT]
        chromosomes = list(range(1, 23)) + ["X", "Y"]
        for bam_filename, output_filename in zip(bam_files, output_files):
            for chromosome in chromosomes:
                command = cls._create_command(
                    library_paths=library_paths,
                    bam_filename=bam_filename,
                    chromosome=chromosome,
                    output_filename=output_filename,
                )
                run(command)
