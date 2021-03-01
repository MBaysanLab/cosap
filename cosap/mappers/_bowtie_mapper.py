import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._mappers import _Mappable, _Mapper


class Bowtie2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(cls, fastq_info: Dict, mapper_config: Dict) -> List:
        flags = mapper_config[MappingKeys.PARAMS]
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
        mapping_config: Dict,
        fastq_info: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
    ) -> List:
        # todo: output filename
        read_group = cls._create_read_group(
            fastq_info=fastq_info, mapping_config=mapping_config
        )
        command = [
            "bowtie2",
            "-p",
            app_config.THREADS,
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
            app_config.THREADS,
            "-bS",
            "-",
            ">",
            mapping_config[MappingKeys.OUTPUT],
        ]
        return command

    @classmethod
    def map(cls, mapping_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        for fastq_info in mapping_config[MappingKeys.DATA]:
            command = cls._create_command(
                mapping_config=mapping_config,
                fastq_info=fastq_info,
                library_paths=library_paths,
                app_config=app_config,
            )
            run(command, cwd=pipeline_config.FASTQ_DIR)
