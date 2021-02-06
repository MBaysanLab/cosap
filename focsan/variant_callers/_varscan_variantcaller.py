import os
from subprocess import run
from pathlib import Path
import glob
from typing import Dict, List, Union

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._variantcallers import _Callable, _VariantCaller

class VarScanVariantCaller(_Callable, _VariantCaller):
    
    @classmethod
    def _create_samtools_mpileup_command(cls, pipeline_config=PipelineConfig, library_paths=LibraryPaths) -> list:
        bam_paths = cls._get_bam_paths(pipeline_config)

        tumor_sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])

        snp_output_name = cls._create_output_filename(pipeline_config, sample_name = "SNP/" + tumor_sample_name)
        indel_output_name = cls._create_output_filename(pipeline_config, sample_name = "INDEL/" + tumor_sample_name)

        command = [
            "samtools",
            "mpileup",
            library_paths.REF_DIR,
            "-q",
            "1",
            "-B",
            bam_paths["germline_bam_path"],
            " ",
            bam_paths["tumor_bam_path"],
            "java",
            "-jar",
            library_paths.VARSCAN,
            "somatic",
            "--output-snp",
            snp_output_name,
            "--output-indel",
            indel_output_name,
            "--mpileup",
            "1",
            "--strand-filter",
            "0"
        ]
        return command

    @classmethod
    def _create_processSomatic_command(cls, pipeline_config=PipelineConfig, library_paths=LibraryPaths, mpileup_object=Union[str, Path]) -> list:
        command = [
            "java",
            "-jar",
            library_paths.VARSCAN,
            "processSomatic",
            mpileup_object,
            "--min-tumor-freq",
            "0.1",
            "--max-tumor-freq",
            "0.05",
            "--p-value",
            "0.07"
        ]
        return command
    
    @classmethod
    def call_variants(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()

        samtools_mpileup = cls._create_samtools_mpileup_command(
            pipeline_config=pipeline_config,
            library_paths=library_paths
        )

        run(samtools_mpileup, cwd=pipeline_config.VCF_OUTPUT_DIR)

        samtools_pileups = glob.glob("*" + pipeline_config.VCF_OUTPUT_DIR + "*vcf*")


        for vcf_file in samtools_pileups:
            processSomatic_command = cls._create_processSomatic_command(
                pipeline_config=pipeline_config,
                library_paths=library_paths,
                mpileup_object=vcf_file
            )
            run(processSomatic_command, cwd=pipeline_config.VCF_OUTPUT_DIR)