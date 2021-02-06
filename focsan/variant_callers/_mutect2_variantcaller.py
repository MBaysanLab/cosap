import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._variantCallers import _Callable, _VariantCaller


class Mutect2VariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls, pipeline_config: PipelineConfig, library_paths: LibraryPaths
    ) -> list:
        bam_paths = cls._get_bam_paths(pipeline_config)

        germline_sample_name = cls._get_sample_name(bam_paths["germline_bam_path"])
        tumor_sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])

        output_name = cls._create_output_filename(
            pipeline_config, sample_name=tumor_sample_name
        )

        command = [
            library_paths.GATK4,
            "Mutect2",
            "-R",
            library_paths.REF_DIR,
            "-I",
            bam_paths["tumor_bam_path"],
            "-tumor",
            tumor_sample_name,
            "-I",
            bam_paths["germline_bam_path"],
            "-normal",
            germline_sample_name,
            "-O",
            output_name,
        ]
        return command

    @classmethod
    def _create_get_snp_variants_command(
        cls, pipeline_config: PipelineConfig, library_paths: LibraryPaths
    ) -> str:
        bam_paths = cls._get_bam_paths(pipeline_config)

        tumor_sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])

        input_name = cls._create_output_filename(
            pipeline_config, sample_name=tumor_sample_name
        )
        output_name = cls._create_output_filename(
            pipeline_config, sample_name=f"SNP_{tumor_sample_name}"
        )

        command = [
            library_paths.GATK4,
            "SeleckVariants",
            "-R",
            library_paths.REF_DIR,
            "-V",
            input_name,
            "--select-type-to-include",
            "SNP",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_indel_variants_command(
        cls, pipeline_config: PipelineConfig, library_paths: LibraryPaths
    ) -> str:
        bam_paths = cls._get_bam_paths(pipeline_config)

        tumor_sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])

        input_name = cls._create_output_filename(
            pipeline_config, sample_name=tumor_sample_name
        )
        output_name = cls._create_output_filename(
            pipeline_config, sample_name=f"SNP_{tumor_sample_name}"
        )

        command = [
            library_paths.GATK4,
            "SeleckVariants",
            "-R",
            library_paths.REF_DIR,
            "-V",
            input_name,
            "--select-type-to-include",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_other_variants_command(
        cls, pipeline_config: PipelineConfig, library_paths: LibraryPaths
    ) -> str:
        bam_paths = cls._get_bam_paths(pipeline_config)

        tumor_sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])

        input_name = cls._create_output_filename(
            pipeline_config, sample_name=tumor_sample_name
        )
        output_name = cls._create_output_filename(
            pipeline_config, sample_name=f"SNP_{tumor_sample_name}"
        )

        command = [
            library_paths.GATK4,
            "SeleckVariants",
            "-R",
            library_paths.REF_DIR,
            "-V",
            input_name,
            "--select-type-to-exclude",
            "SNP",
            "--select-type-to-exclude",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def call_variants(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()

        mutect_command = cls._create_run_command(
            pipeline_config=pipeline_config, library_paths=library_paths
        )
        get_snp_command = cls._create_get_snp_variants_command(
            pipeline_config=pipeline_config, library_paths=library_paths
        )
        get_indel_command = cls._create_get_indel_variants_command(
            pipeline_config=pipeline_config, library_paths=library_paths
        )
        get_other_variants_command = cls._create_get_other_variants_command(
            pipeline_config=pipeline_config, library_paths=library_paths
        )

        run(mutect_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
        run(get_snp_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
        run(get_indel_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
        run(get_other_variants_command, cwd=pipeline_config.VCF_OUTPUT_DIR)
