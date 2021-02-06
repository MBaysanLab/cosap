import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._variantCallers import _Callable, _VariantCaller


class Strelka2VariantCaller(_Callable, _VariantCaller):
    @classmethod
<<<<<<< HEAD
    def _create_strelka_command(
        cls, pipeline_config=PipelineConfig, library_paths=LibraryPaths
    ) -> List:
=======
    def _create_strelka_command(cls, pipeline_config=PipelineConfig, library_paths=LibraryPaths) -> List:
>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a
        bam_paths = cls._get_bam_paths(pipeline_config)
        command = [
            library_paths.STRELKA,
            "--normalBam",
            bam_paths["germline_bam_path"],
            "--tumorBam",
            bam_paths["tumor_bam_path"],
            "--referenceFasta",
            library_paths.REF_DIR,
            "--runDir",
            pipeline_config.VCF_OUTPUT_DIR,
            "--exome",
<<<<<<< HEAD
            "--disableEVS",
        ]

        return command

=======
            "--disableEVS"
        ]

        return command
    
>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a
    @classmethod
    def call_variants(cls, pipeline_config=PipelineConfig):
        library_paths = LibraryPaths()
        strelka_command = cls._create_strelka_command(
<<<<<<< HEAD
            pipeline_config=pipeline_config, library_paths=library_paths
        )

        run(strelka_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
=======
            pipeline_config=pipeline_config,
            library_paths=library_paths
        )

        run(strelka_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
    
    
>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a
