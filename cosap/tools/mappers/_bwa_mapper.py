from subprocess import PIPE, Popen, check_output
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MappingKeys
from ._mappers import _Mappable, _Mapper


class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(cls, mapper_config: Dict) -> str:
        flags = cls._create_readgroup_flags(
            mapper_config=mapper_config,
        )

        read_arguments = []
        if MappingKeys.RG_ID in flags.keys():
            read_arguments.append(fr"@RG\tID:{flags[MappingKeys.RG_ID]}")
        if MappingKeys.RG_SM in flags.keys():
            read_arguments.append(fr"@RG\tSM:{flags[MappingKeys.RG_SM]}")
        if MappingKeys.RG_LB in flags.keys():
            read_arguments.append(fr"@RG\tLB:{flags[MappingKeys.RG_LB]}")
        if MappingKeys.RG_PL in flags.keys():
            read_arguments.append(fr"@RG\tPL:{flags[MappingKeys.RG_PL]}")
        if MappingKeys.RG_PU in flags.keys():
            read_arguments.append(fr"@RG\tPU:{flags[MappingKeys.RG_PU]}")

        read_groups = "".join(read_arguments)
        return read_groups

    @classmethod
    def _create_command(
        cls,
        mapper_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        read_group: str = None,
    ) -> List:
        fastq_inputs = [fastq for fastq in mapper_config[MappingKeys.INPUT].values()]

        command = [
            "bwa",
            "mem",
            "-t",
            str(app_config.MAX_THREADS_PER_JOB),
            library_paths.BWA_ASSEMBLY,
            *fastq_inputs,
        ]

        if read_group:
            command.extend(["-R", read_group])
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
        sort_command = cls._samtools_sort_command(
            app_config=app_config, output_path=mapper_config[MappingKeys.OUTPUT]
        )
        index_command = cls._samtools_index_command(
            app_config=app_config, input_path=mapper_config[MappingKeys.OUTPUT]
        )

        bwa = Popen(bwa_command, stdout=PIPE)
        samtools = check_output(sort_command, stdin=bwa.stdout)
        bwa.wait()
        # run(index_command)
