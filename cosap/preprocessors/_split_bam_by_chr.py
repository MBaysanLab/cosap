from subprocess import run
from typing import Dict, List
from pathlib import Path

from .._library_paths import LibraryPaths
from .._pipeline_config import SplitKeys

class SplitbyCHR:
    
    @classmethod
    def _create_command(cls, library_paths:LibraryPaths, bam_file:Path, chr:str) -> List:
        
        command = [
            "samtools",
            "view",
            "-bh",
            bam_file,
            f"chr{chr}",
            ">",
            f"{bam_file}_CHR_{chr}.bam"
        ]
        return command
    
    @classmethod
    def split_by_chr(cls, library_paths:LibraryPaths, split_config:Dict):
        bam_files = split_config[SplitKeys.INPUT]
        chromosomes = [i for i in range(1,23)] + ["X", "Y"]
        for bam in bam_files:
            for chrom in chromosomes:
                command = cls._create_command(
                                library_paths=library_paths,
                                bam_file=bam,
                                chr=chrom
                            )
                run(command, cwd=split_config.OUTPUT_DIR)
        