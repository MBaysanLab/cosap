import os
from pathlib import Path
from subprocess import PIPE, Popen, check_output
from typing import Dict, List

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MappingKeys
from ...pipeline_runner.runners import DockerRunner
from ._mappers import _Mappable, _Mapper


class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, mapper_config: Dict, is_parabricks: bool = False
    ) -> str:
        flags = cls._create_readgroup_flags(
            mapper_config=mapper_config,
        )

        read_arguments = []

        if not is_parabricks:
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

        else:
            if MappingKeys.RG_ID in flags.keys():
                read_arguments.append(
                    rf"--read-group-id-prefix {flags[MappingKeys.RG_ID]}"
                )
            if MappingKeys.RG_SM in flags.keys():
                read_arguments.append(rf"--read-group-sm {flags[MappingKeys.RG_SM]}")
            if MappingKeys.RG_LB in flags.keys():
                read_arguments.append(rf"--read-group-lb {flags[MappingKeys.RG_LB]}")
            if MappingKeys.RG_PL in flags.keys():
                read_arguments.append(rf"--read-group-pl {flags[MappingKeys.RG_PL]}")

            read_groups = " ".join(read_arguments)

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
    def _create_parabricks_command(
        cls,
        mapper_config: Dict,
        library_paths: LibraryPaths,
        read_group: str = None,
    ) -> List:

        fastq_inputs = [fastq for fastq in mapper_config[MappingKeys.INPUT].values()]

        command = [
            "pbrun",
            "fq2bam",
            "--ref",
            library_paths.BWA_ASSEMBLY,
            "--in-fq",
            *fastq_inputs,
            "--out-bam",
            mapper_config[MappingKeys.OUTPUT],
            "--no-markdups",
        ]

        if read_group:
            command.extend([read_group])
        return command

    @classmethod
    def map(cls, mapper_config: Dict, device: str = "cpu") -> None:
        library_paths = LibraryPaths()
        app_config = AppConfig()
        if device == "cpu":
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

        elif device == "gpu":
            read_group = cls._create_read_group(
                mapper_config=mapper_config, is_parabricks=True
            )

            bwa_command = cls._create_parabricks_command(
                mapper_config=mapper_config,
                read_group=read_group,
                library_paths=library_paths,
            )
            output_dir = os.path.abspath(
                os.path.dirname(mapper_config[MappingKeys.OUTPUT])
            )
            os.makedirs(output_dir, exist_ok=True)
            runner = DockerRunner(device=device)
            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(bwa_command),
                workdir=str(Path(output_dir).parent.parent),
            )
