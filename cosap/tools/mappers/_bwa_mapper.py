import os
from pathlib import Path
from subprocess import PIPE, Popen
from typing import Dict, List, Union

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MappingKeys
from ..._utils import convert_to_absolute_path
from ...runners.docker_runner import DockerRunner
from ._mappers import _Mappable, _Mapper


class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, mapper_config: Dict, is_parabricks: bool = False
    ) -> Union[str, List[str]]:
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
                read_arguments.append("--read-group-id-prefix")
                read_arguments.append(str(flags[MappingKeys.RG_ID]))
            if MappingKeys.RG_SM in flags.keys():
                read_arguments.extend(
                    ["--read-group-sm", str(flags[MappingKeys.RG_SM])]
                )
            if MappingKeys.RG_LB in flags.keys():
                read_arguments.extend(
                    ["--read-group-lb", str(flags[MappingKeys.RG_LB])]
                )
            if MappingKeys.RG_PL in flags.keys():
                read_arguments.extend(
                    ["--read-group-pl", str(flags[MappingKeys.RG_PL])]
                )

            read_groups = read_arguments

        return read_groups

    @classmethod
    def _create_command(
        cls,
        mapper_config: Dict,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        read_group: str = None,
    ) -> List:
        fastq_inputs = [
            convert_to_absolute_path(fastq)
            for fastq in mapper_config[MappingKeys.INPUT].values()
        ]

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

    @staticmethod
    def _create_parabricks_output_path(mapper_config: Dict) -> str:
        output_path = mapper_config[MappingKeys.OUTPUT]
        if os.path.isabs(output_path):
            return output_path
        return os.path.join(mapper_config[MappingKeys.OUTPUT_DIR], output_path)

    @classmethod
    def _create_parabricks_command(
        cls,
        mapper_config: Dict,
        library_paths: LibraryPaths,
        read_group: List[str] = None,
    ) -> List:
        fastq_inputs = [
            convert_to_absolute_path(fastq)
            for fastq in mapper_config[MappingKeys.INPUT].values()
        ]
        output_path = cls._create_parabricks_output_path(mapper_config)

        command = [
            "pbrun",
            "fq2bam",
            "--ref",
            library_paths.BWA_ASSEMBLY,
            "--in-fq",
            *fastq_inputs,
            "--out-bam",
            output_path,
            "--bwa-options=-K 10000000",
        ]
        if not mapper_config[MappingKeys.POST_PROCESSING]:
            command.extend(["--no-markdups"])
        else:
            command.extend(
                [
                    "--knownSites",
                    library_paths.DBSNP,
                    "--knownSites",
                    library_paths.MILLS_INDEL,
                    "--knownSites",
                    library_paths.ONE_THOUSAND_G,
                    "--out-recal-file",
                    f"{output_path}.recal",
                ]
            )

        if read_group:
            command.extend(read_group)
        return command

    @classmethod
    def map(cls, mapper_config: Dict, device: str = "cpu") -> None:
        library_paths = LibraryPaths()
        app_config = AppConfig()
        workdir = mapper_config[MappingKeys.OUTPUT_DIR]
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
            # index_command = cls._samtools_index_command(
            #     app_config=app_config, input_path=mapper_config[MappingKeys.OUTPUT]
            # )

            bwa = Popen(bwa_command, stdout=PIPE, stderr=PIPE, cwd=workdir)
            samtools = Popen(
                sort_command,
                stdin=bwa.stdout,
                stdout=PIPE,
                stderr=PIPE,
                cwd=workdir,
            )

            bwa.stdout.close()
            _, samtools_stderr = samtools.communicate()
            _, bwa_stderr = bwa.communicate()

            if bwa.returncode != 0:
                raise Exception(
                    f"BWA command {(*bwa_command,)} failed with error: {bwa_stderr.decode('utf-8')}"
                )
            if samtools.returncode != 0:
                raise Exception(
                    f"samtools sort command {(*sort_command,)} failed with error: {samtools_stderr.decode('utf-8')}"
                )
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
            input_dir = os.path.abspath(
                os.path.dirname(mapper_config[MappingKeys.INPUT]["1"])
            )
            output_path = cls._create_parabricks_output_path(mapper_config)
            output_dir = os.path.dirname(output_path)
            if not os.path.isabs(output_dir):
                output_dir = os.path.join(
                    str(Path(mapper_config[MappingKeys.WORKDIR]).absolute()), output_dir
                )
            os.makedirs(output_dir, exist_ok=True)
            runner = DockerRunner(device=device)
            runner.run(
                image=DockerImages.PARABRICKS,
                command=bwa_command,
                workdir=str(Path(mapper_config[MappingKeys.WORKDIR]).absolute()),
                paths_to_bind=[
                    {"source": input_dir, "target": input_dir},
                ],
            )
