import os
from pathlib import Path
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MDUPKeys
from ..._utils import convert_to_absolute_path
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
        tmp_dir = memory_handler.get_temp_dir(
            dir=os.path.dirname(mdup_config[MDUPKeys.OUTPUT_DIR])
        )

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
        tmp_dir = memory_handler.get_temp_dir(
            dir=os.path.dirname(mdup_config[MDUPKeys.OUTPUT_DIR])
        )

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
        input_bam: str = None,
        output_bam: str = None,
        metrics_path: str = None,
    ) -> List:
        input_bam = input_bam or mdup_config[MDUPKeys.INPUT]
        output_bam = output_bam or cls._create_parabricks_output_path(mdup_config)
        metrics_path = metrics_path or f"{output_bam}.duplicate_metrics"

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
            metrics_path,
        ]

        return command

    @staticmethod
    def _create_parabricks_output_path(mdup_config: Dict) -> str:
        output_path = mdup_config[MDUPKeys.OUTPUT]
        if os.path.isabs(output_path):
            return output_path
        return os.path.join(mdup_config[MDUPKeys.OUTPUT_DIR], output_path)

    @classmethod
    def _create_parabricks_intermediate_path(
        cls, mdup_config: Dict, suffix: str
    ) -> str:
        output_path = cls._create_parabricks_output_path(mdup_config)
        root, extension = os.path.splitext(output_path)
        return f"{root}.{suffix}{extension}"

    @classmethod
    def create_parabricks_sort_command(
        cls,
        input_bam: str,
        output_bam: str,
        library_paths: LibraryPaths,
        sort_order: str,
    ) -> list:
        command = [
            "pbrun",
            "bamsort",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            input_bam,
            "--out-bam",
            output_bam,
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
                run(command, cwd=mdup_config[MDUPKeys.OUTPUT_DIR], check=True)

        elif device == "gpu":
            output_bam = cls._create_parabricks_output_path(mdup_config)
            queryname_bam = cls._create_parabricks_intermediate_path(
                mdup_config, "queryname"
            )
            markdup_bam = cls._create_parabricks_intermediate_path(mdup_config, "markdup")

            output_dir = os.path.dirname(output_bam)
            if not os.path.isabs(output_dir):
                output_dir = os.path.join(
                    str(Path(mdup_config[MDUPKeys.WORKDIR]).absolute()), output_dir
                )
            os.makedirs(output_dir, exist_ok=True)

            runner = DockerRunner(device=device)
            # Parabricks markdup requires input bam to be queryname sorted
            query_name_sort_command = cls.create_parabricks_sort_command(
                mdup_config[MDUPKeys.INPUT], queryname_bam, library_paths, "queryname"
            )
            return_code = runner.run(
                image=DockerImages.PARABRICKS,
                command=query_name_sort_command,
                workdir=str(Path(mdup_config[MDUPKeys.WORKDIR]).absolute()),
            )
            if return_code not in (0, None):
                raise RuntimeError(
                    f"Parabricks queryname sort failed with return code {return_code}"
                )

            mdup_command = cls.create_parabricks_markdup_command(
                library_paths=library_paths,
                app_config=app_config,
                mdup_config=mdup_config,
                memory_handler=None,
                input_bam=queryname_bam,
                output_bam=markdup_bam,
                metrics_path=f"{output_bam}.duplicate_metrics",
            )

            return_code = runner.run(
                image=DockerImages.PARABRICKS,
                command=mdup_command,
                workdir=str(Path(mdup_config[MDUPKeys.WORKDIR]).absolute()),
            )
            if return_code not in (0, None):
                raise RuntimeError(
                    f"Parabricks markdup failed with return code {return_code}"
                )

            # Convert output to coordinate sorted bam
            coordinate_sort_command = cls.create_parabricks_sort_command(
                markdup_bam, output_bam, library_paths, "coordinate"
            )
            return_code = runner.run(
                image=DockerImages.PARABRICKS,
                command=coordinate_sort_command,
                workdir=str(Path(mdup_config[MDUPKeys.WORKDIR]).absolute()),
            )
            if return_code not in (0, None):
                raise RuntimeError(
                    f"Parabricks coordinate sort failed with return code {return_code}"
                )
