import os
from posixpath import commonpath
from subprocess import PIPE, STDOUT, Popen, check_output, run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MappingKeys
from ._mappers import _Mappable, _Mapper


class BWA2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(cls, mapper_config: Dict) -> str:
        flags = mapper_config[MappingKeys.PARAMS][MappingKeys.READ_GROUP]
        read_arguments = "".join(
            (
                r"@RG\tID:",
                flags[MappingKeys.RG_ID],
                r"\tSM:",
                flags[MappingKeys.RG_SM],
                r"\tLB:",
                flags[MappingKeys.RG_LB],
                r"\tPL:",
                flags[MappingKeys.RG_PL],
                r"\tPU:",
                flags[MappingKeys.RG_PU],
            )
        )
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        mapper_config: Dict,
        read_group: str,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:

        fastq_inputs = [fastq for fastq in mapper_config[MappingKeys.INPUT].values()]

        command = [
            "bwa-mem2",
            "mem",
            "-t",
            str(app_config.THREADS),
            "-R",
            read_group,
            library_paths.BWA_ASSEMBLY,
            *fastq_inputs,
        ]
        return command

    @classmethod
    def _create_samtools_command(
        cls,
        mapper_config: Dict,
        read_group: str,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:

        command = [
            "samtools",
            "sort",
            "-@",
            str(app_config.THREADS),
            "-o",
            mapper_config[MappingKeys.OUTPUT],
            "-",
        ]
        return command

    @classmethod
    def map(cls, mapper_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        read_group = cls._create_read_group(mapper_config=mapper_config)

        bwa_command = cls._create_command(
            mapper_config=mapper_config,
            read_group=read_group,
            library_paths=library_paths,
            app_config=app_config,
        )
        samtools_command = cls._create_samtools_command(
            mapper_config=mapper_config,
            read_group=read_group,
            library_paths=library_paths,
            app_config=app_config,
        )
        bwa = Popen(bwa_command, stdout=PIPE)
        samtools = check_output(samtools_command, stdin=bwa.stdout)
        bwa.wait()