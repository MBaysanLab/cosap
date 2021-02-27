import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._mappers import _Mappable, _Mapper


class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, fastq_info: Dict, mapper_config: Dict
    ) -> str:
        flags = mapper_config[MappingKeys.PARAMS]
        read_arguments = "".join(
            (
                '"',
                r"@RG\tID",
                flags[MappingKeys.RG_ID],
                r"\tSM:",
                flags[MappingKeys.RG_SM],
                r"\tLB",
                flags[MappingKeys.RG_LB],
                r"\tPL",
                flags[MappingKeys.RG_PL],
                r"\tPU",
                flags[MappingKeys.RG_PU],
                '"',
            )
        )
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        mapper_config: Dict,
        fastq_info: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:
        output_filename = cls._create_output_filename(
            fastq_info=fastq_info, mapper_config=mapper_config
        )
        read_group = cls._create_read_group(
            fastq_info=fastq_info, mapper_config=mapper_config
        )
        command = [
            "bwa",
            "mem",
            "-t",
            app_config.THREADS,
            "-R",
            read_group,
            library_paths.REF_DIR,
            # TODO: This needs to go somewhere else
            os.path.normpath("Bwa/Homo_sapiens_assembly38.fasta"),
            fastq_info["R1"],
            fastq_info["R2"],
            "|",
            "samtools",
            "view",
            "-@",
            app_config.THREADS,
            "-bS",
            "-",
            ">",
            output_filename,
        ]
        return command

    @classmethod
    def map(cls, mapper_config: PipelineConfig):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        for fastq_info in mapping_config[MappingKeys.DATA]:
            command = cls._create_command(
                mapper_config=mapper_config,
                fastq_info=fastq_info,
                library_paths=library_paths,
                app_config=app_config,
            )
            run(command, cwd=mapper_config.FASTQ_DIR)
