import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._variantCallers import _Callable, _VariantCaller

<<<<<<< HEAD

class SomaticSniperVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_somaticSniper_command(
        cls, pipeline_config=PipelineConfig, library_paths=LibraryPaths
    ) -> List:
        bam_paths = cls._get_bam_paths(pipeline_config)

        sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])
        output_name = cls._create_output_filename(
            pipeline_config, sample_name=sample_name
        )
=======
class SomaticSniperVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_somaticSniper_command(cls, pipeline_config=PipelineConfig, library_paths = LibraryPaths) -> List:
        bam_paths = cls._get_bam_paths(pipeline_config)

        sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])
        output_name = cls._create_output_filename(pipeline_config, sample_name=sample_name)
>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a

        command = [
            library_paths.SOMATICSNIPER,
            "-q",
            "1",
            "-L",
            "-G",
            "-Q",
            "15",
            "-s",
            "0.01",
            "-T",
            "0.85",
            "-N",
            "2",
            "-r",
            "0.001",
            "-n",
            "NORMAL",
            "-t",
            "TUMOR",
            "-F",
            "vcf",
            "-f",
            library_paths.REF_DIR,
            bam_paths["tumor_bam_path"],
            bam_paths["germline_bam_path"],
<<<<<<< HEAD
            output_name,
=======
            output_name
>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a
        ]

        return command

<<<<<<< HEAD
=======

>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a
    @classmethod
    def call_variants(cls, pipeline_config=PipelineConfig):
        library_paths = LibraryPaths()

        somatic_sniper_command = cls._create_somaticSniper_command(
<<<<<<< HEAD
            pipeline_config=pipeline_config, library_paths=library_paths
        )
        run(somatic_sniper_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
=======
            pipeline_config=pipeline_config,
            library_paths=library_paths
        )
        run(somatic_sniper_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
>>>>>>> 55f6e9fd8145ae365b07eea60a3395f9cb216f5a
