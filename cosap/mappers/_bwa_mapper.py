import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MappingKeys
from ._mappers import _Mappable, _Mapper


class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(cls, mapper_config: Dict) -> str:
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
        read_group: str,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:
        command = [
            "bwa",
            "mem",
            "-t",
            app_config.THREADS,
            "-R",
            read_group,
            library_paths.REF_DIR,
            # TODO: This needs to go somewhere else
            library_paths.BWA_ASSEMBLY,
            *mapper_config[MappingKeys.INPUTS],
            "|",
            "samtools",
            "view",
            "-@",
            app_config.THREADS,
            "-bS",
            "-",
            ">",
            mapper_config[MappingKeys.OUTPUT],
        ]
        return command

    @classmethod
    def map(cls, mapper_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        read_group = cls._create_read_group(mapper_config=mapper_config)

        command = cls._create_command(
            mapper_config=mapper_config,
            read_group=read_group,
            library_paths=library_paths,
            app_config=app_config,
        )

        run(command, cwd=os.path.dirname(mapper_config[MappingKeys.INPUTS][0]))
