import os
from subprocess import Popen, run, check_output, PIPE
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MappingKeys
from ._mappers import _Mappable, _Mapper


class Bowtie2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_fastq_reads_command(cls, mapper_config: Dict) -> List:
        command = []
        for i, read in enumerate(mapper_config[MappingKeys.INPUT], 1):
            command.extend([f"-{i}", mapper_config[MappingKeys.INPUT][read]])
        return command

    @classmethod
    def _create_read_group(cls, mapper_config: Dict) -> List:
        flags = mapper_config[MappingKeys.PARAMS][MappingKeys.READ_GROUP]
        read_arguments = [
            "--rg-id",
            flags[MappingKeys.RG_ID],
            "--rg",
            f"SM:{flags[MappingKeys.RG_SM]}",
            "--rg",
            f"LB:{flags[MappingKeys.RG_LB]}",
            "--rg",
            f"PL:{flags[MappingKeys.RG_PL]}",
            "--rg",
            f"PU:{flags[MappingKeys.RG_PU]}",
        ]
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        mapper_config: Dict,
        read_group: List,
        fastq_reads: List,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:
        command = [
            "bowtie2",
            "-p",
            str(app_config.THREADS),
            *read_group,
            "-x",
            library_paths.BOWTIE2_ASSEMBLY,
            *fastq_reads,
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

        fastq_reads = cls._create_fastq_reads_command(mapper_config=mapper_config)

        bowtie_command = cls._create_command(
            mapper_config=mapper_config,
            read_group=read_group,
            fastq_reads=fastq_reads,
            library_paths=library_paths,
            app_config=app_config,
        )

        samtools_command = cls._create_samtools_command(
            mapper_config=mapper_config,
            read_group=read_group,
            library_paths=library_paths,
            app_config=app_config,
        )

        bowtie = Popen(bowtie_command, stdout=PIPE)
        samtools = check_output(samtools_command, stdin=bowtie.stdout)
        bowtie.wait()
