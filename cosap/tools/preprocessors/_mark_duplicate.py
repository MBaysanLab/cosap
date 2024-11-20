import os
from pathlib import Path
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MDUPKeys
from ..._utils import convert_to_absolute_path, join_paths
from ...memory_handler import MemoryHandler
from ...runners.docker_runner import DockerRunner
from ._preprocessors import _PreProcessable, _Preprocessor


class MarkDuplicate(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command_spark(
        cls,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        mdup_config: Dict,
        memory_handler: MemoryHandler,
    ) -> List:
        input_bam = convert_to_absolute_path(
            memory_handler.get_path(mdup_config[MDUPKeys.INPUT])
        )
        tmp_dir = memory_handler.get_temp_dir(dir=os.path.dirname(MDUPKeys.OUTPUT_DIR))

        command = [
            "gatk",
            "MarkDuplicatesSpark",
            "-I",
            input_bam,
            "-O",
            mdup_config[MDUPKeys.OUTPUT],
            "-M",
            f"{mdup_config[MDUPKeys.OUTPUT]}_metrics",
            "--create-output-bam-index",
            "--spark-master",
            f"local[{app_config.MAX_THREADS_PER_JOB}]",
            # "--tmp-dir",
            # tmp_dir,
            "--verbosity",
            "WARNING",
        ]
        if mdup_config[MDUPKeys.DUPLICATE_HANDLING_METHOD] == "delete":
            command.append("--remove-all-duplicates")

        return command

    @classmethod
    def _create_command(
        cls,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        mdup_config: Dict,
        memory_handler: MemoryHandler,
    ) -> List:
        input_bam = convert_to_absolute_path(
            memory_handler.get_path(mdup_config[MDUPKeys.INPUT])
        )
        tmp_dir = memory_handler.get_temp_dir(dir=os.path.dirname(MDUPKeys.OUTPUT_DIR))

        command = [
            "picard",
            "MarkDuplicates",
            "--INPUT",
            input_bam,
            "--OUTPUT",
            mdup_config[MDUPKeys.OUTPUT],
            "--METRICS_FILE",
            f"{mdup_config[MDUPKeys.OUTPUT]}_metrics",
            "--CREATE_INDEX",
            "true",
            "--TMP_DIR",
            tmp_dir,
        ]
        if mdup_config[MDUPKeys.DUPLICATE_HANDLING_METHOD] == "delete":
            command.append("--REMOVE_DUPLICATES")
            command.append("true")
        return command

    @classmethod
    def create_parabricks_markdup_command(
        cls,
        library_paths: LibraryPaths,
        app_config: AppConfig,
        mdup_config: Dict,
        memory_handler: MemoryHandler,
    ) -> List:
        input_bam = mdup_config[MDUPKeys.INPUT]
        output_bam = mdup_config[MDUPKeys.OUTPUT]

        command = [
            "pbrun",
            "markdup",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            input_bam,
            "--out-bam",
            output_bam,
            "--out-duplicate-metrics",
            f"{output_bam}.duplicate_metrics",
        ]

        return command

    @classmethod
    def create_parabricks_sort_command(
        cls, mdup_config: dict, library_paths: LibraryPaths, sort_order: str
    ) -> list:

        command = [
            "pbrun",
            "bamsort",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            mdup_config[MDUPKeys.INPUT],
            "--out-bam",
            mdup_config[MDUPKeys.INPUT],
            "--sort-order",
            sort_order,
        ]

        return command

    @classmethod
    def run_preprocessor(cls, mdup_config: Dict, device: str = "cpu"):
        app_config = AppConfig()
        library_paths = LibraryPaths()

        if device == "cpu":
            command_func = (
                cls._create_command_spark
                if mdup_config[MDUPKeys.SPARK]
                else cls._create_command
            )

            with MemoryHandler() as memory_handler:
                command = command_func(
                    library_paths=library_paths,
                    app_config=app_config,
                    mdup_config=mdup_config,
                    memory_handler=memory_handler,
                )
                run(command, cwd=mdup_config[MDUPKeys.OUTPUT_DIR])

        elif device == "gpu":

            output_dir = os.path.abspath(os.path.dirname(mdup_config[MDUPKeys.OUTPUT]))
            os.makedirs(output_dir, exist_ok=True)

            runner = DockerRunner(device=device)
            # Parabricks markdup requires input bam to be queryname sorted
            query_name_sort_command = cls.create_parabricks_sort_command(
                mdup_config, library_paths, "queryname"
            )
            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(query_name_sort_command),
                workdir=str(Path(output_dir).parent.parent),
            )

            mdup_command = cls.create_parabricks_markdup_command(
                library_paths=library_paths,
                app_config=app_config,
                mdup_config=mdup_config,
                memory_handler=None,
            )

            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(mdup_command),
                workdir=str(Path(output_dir).parent.parent),
            )

            # Convert output to coordinate sorted bam
            coordinate_sort_command = cls.create_parabricks_sort_command(
                mdup_config, library_paths, "coordinate"
            )
            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(coordinate_sort_command),
                workdir=str(Path(output_dir).parent.parent),
            )
