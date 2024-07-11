from subprocess import PIPE, Popen, check_output
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MappingKeys
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
        flags = cls._create_readgroup_flags(
            mapper_config=mapper_config,
        )
        read_arguments = []
        if MappingKeys.RG_ID in flags.keys():
            read_arguments.append(rf"--rg-id {flags[MappingKeys.RG_ID]}")
        if MappingKeys.RG_SM in flags.keys():
            read_arguments.append(rf"--rg SM:{flags[MappingKeys.RG_SM]}")
        if MappingKeys.RG_LB in flags.keys():
            read_arguments.append(rf"--rg LB:{flags[MappingKeys.RG_LB]}")
        if MappingKeys.RG_PL in flags.keys():
            read_arguments.append(rf"--rg PL:{flags[MappingKeys.RG_PL]}")
        if MappingKeys.RG_PU in flags.keys():
            read_arguments.append(rf"--rg PU:{flags[MappingKeys.RG_PU]}")

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
            str(app_config.MAX_THREADS_PER_JOB),
            "-x",
            library_paths.BOWTIE2_ASSEMBLY,
            *fastq_reads,
            *read_group,
        ]

        return command

    @classmethod
    def map(cls, mapper_config: Dict, *args, **kwargs):
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

        sort_command = cls._samtools_sort_command(
            app_config=app_config, output_path=mapper_config[MappingKeys.OUTPUT]
        )

        index_command = cls._samtools_index_command(
            app_config=app_config, input_path=mapper_config[MappingKeys.OUTPUT]
        )

        bowtie = Popen(bowtie_command, stdout=PIPE, cwd=mapper_config[MappingKeys.OUTPUT_DIR])
        samtools = check_output(sort_command, stdin=bowtie.stdout, cwd=mapper_config[MappingKeys.OUTPUT_DIR])
        bowtie.wait()

        if bowtie.returncode != 0:
            raise Exception("Bowtie2 mapper failed.")
        else:
            print(samtools.decode("utf-8"))
        # run(index_command)
