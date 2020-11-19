import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._mappers import _Mappable, _Mapper


class Bowtie2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, fastq_info: Dict, pipeline_config: PipelineConfig
    ) -> List:
        flags = cls._get_rg_flags(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_arguments = [
            "--rg-id",
            flags[cls.RG_ID],
            "--rg",
            f"SM:{flags[cls.RG_SM]}",
            "--rg",
            f"LB:{flags[cls.RG_LB]}",
            "--rg",
            f"PL:{flags[cls.RG_PL]}",
            "--rg",
            f"PU:{flags[cls.RG_PU]}",
        ]
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        pipeline_config: PipelineConfig,
        fastq_info: Dict,
        library_paths: LibraryPaths,
    ) -> List:
        output_filename = cls._create_output_filename(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_group = cls._create_read_group(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        command = [
            "bowtie2",
            "-p",
            pipeline_config.MAPPER_THREADS,
            *read_group,
            "-x",
            library_paths.REF_DIR,
            # TODO: This needs to go somewhere else
            os.path.normpath("Bowtie2/Homo_sapiens_assembly38"),
            "-1",
            fastq_info["R1"],
            "-2",
            fastq_info["R2"],
            "|",
            "samtools",
            "view",
            "-@",
            pipeline_config.MAPPER_THREADS,
            "-bS",
            "-",
            ">",
            output_filename,
        ]
        return command

    @classmethod
    def map(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()
        fastq_info_list = cls._get_file_information(pipeline_config=pipeline_config)

        for fastq_info in fastq_info_list:
            command = cls._create_command(
                pipeline_config=pipeline_config,
                fastq_info=fastq_info,
                library_paths=library_paths,
            )
            run(command, cwd=pipeline_config.FASTQ_DIR)
