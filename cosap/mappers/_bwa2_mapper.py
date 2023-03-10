from subprocess import PIPE, STDOUT, Popen, check_output, run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MappingKeys
from ._mappers import _Mappable, _Mapper


class BWA2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(cls, mapper_config: Dict) -> str:

        if not MappingKeys.READ_GROUP in mapper_config[MappingKeys.PARAMS].keys():
            return ""

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
            str(app_config.MAX_THREADS_PER_JOB),
            library_paths.BWA_ASSEMBLY,
            *fastq_inputs,
        ]
        if MappingKeys.READ_GROUP in mapper_config[MappingKeys.PARAMS].keys():
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
        #run(index_command)
