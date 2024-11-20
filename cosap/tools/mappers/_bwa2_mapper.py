from subprocess import PIPE, run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MappingKeys
from ..._utils import convert_to_absolute_path
from ._mappers import _Mappable, _Mapper


class BWA2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(cls, mapper_config: Dict) -> str:
        flags = cls._create_readgroup_flags(
            mapper_config=mapper_config,
        )

        read_arguments = []
        if MappingKeys.RG_ID in flags.keys():
            read_arguments.append(rf"\tID:{flags[MappingKeys.RG_ID]}")
        if MappingKeys.RG_SM in flags.keys():
            read_arguments.append(rf"\tSM:{flags[MappingKeys.RG_SM]}")
        if MappingKeys.RG_LB in flags.keys():
            read_arguments.append(rf"\tLB:{flags[MappingKeys.RG_LB]}")
        if MappingKeys.RG_PL in flags.keys():
            read_arguments.append(rf"\tPL:{flags[MappingKeys.RG_PL]}")
        if MappingKeys.RG_PU in flags.keys():
            read_arguments.append(rf"\tPU:{flags[MappingKeys.RG_PU]}")

        if len(read_arguments) > 0:
            read_arguments.insert(0, "@RG")

        read_groups = "".join(read_arguments)

        return read_groups

    @classmethod
    def _create_command(
        cls,
        mapper_config: Dict,
        read_group: str,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:
        fastq_inputs = [
            convert_to_absolute_path(fastq)
            for fastq in mapper_config[MappingKeys.INPUT].values()
        ]

        command = [
            "bwa-mem2",
            "mem",
            "-t",
            str(app_config.MAX_THREADS_PER_JOB),
            library_paths.BWA_ASSEMBLY,
            *fastq_inputs,
        ]
        if MappingKeys.READ_GROUP in mapper_config[MappingKeys.PARAMS].keys():
            command.extend(["-R", read_group])

        return command

    @classmethod
    def map(cls, mapper_config: Dict, *args, **kwargs):
        library_paths = LibraryPaths()
        app_config = AppConfig()
        workdir = mapper_config[MappingKeys.OUTPUT_DIR]
        read_group = cls._create_read_group(mapper_config=mapper_config)

        bwa_command = cls._create_command(
            mapper_config=mapper_config,
            read_group=read_group,
            library_paths=library_paths,
            app_config=app_config,
        )
        sort_command = cls._samtools_sort_command(
            app_config=app_config, output_path=mapper_config[MappingKeys.OUTPUT]
        )
        index_command = cls._samtools_index_command(
            app_config=app_config, input_path=mapper_config[MappingKeys.OUTPUT]
        )
        bwa = run(bwa_command, stdout=PIPE, cwd=workdir)
        samtools = run(
            sort_command,
            input=bwa.stdout,
            stdout=PIPE,
            stderr=PIPE,
            cwd=workdir,
            check=True,
        )

        if bwa.returncode != 0:
            raise Exception(f"BWA2 command {*bwa_command,} failed.")
        # run(index_command)
