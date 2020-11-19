import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._mappers import _Mappable, _Mapper


class NovoalignMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, fastq_info: Dict, pipeline_config: PipelineConfig
    ) -> str:
        flags = cls._get_rg_flags(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_arguments = "".join(
            (
                '"',
                r"@RG\tID",
                flags[cls.RG_ID],
                r"\tSM:",
                flags[cls.RG_SM],
                r"\tLB",
                flags[cls.RG_LB],
                r"\tPL",
                flags[cls.RG_PL],
                r"\tPU",
                flags[cls.RG_PU],
                '"',
            )
        )
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
        stats_filename = f"{os.path.splitext(output_filename)[0]}_stats.txt"
        command = [
            "novoalign",
            "-k",
            "-d",
            library_paths.REF_DIR,
            os.path.normpath("Bowtie2/Homo_sapiens_assembly38"),
            "-f",
            fastq_info["R1"],
            fastq_info["R2"],
            "-a",
            "-c",
            pipeline_config.MAPPER_THREADS,
            "-o",
            "SAM",
            *read_group,
            "2>",
            stats_filename,
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
