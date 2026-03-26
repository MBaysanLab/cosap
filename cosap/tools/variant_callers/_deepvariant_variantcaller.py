import os
from typing import Dict

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ..._utils import convert_to_absolute_path, join_paths
from ...runners.docker_runner import DockerRunner
from ._variantcallers import _Callable, _VariantCaller


class DeepVariantVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def create_run_deepvariant_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        output_vcf = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT],
        )
        output_gvcf = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            caller_config[VariantCallingKeys.GVCF_OUTPUT],
        )

        command = [
            "run_deepvariant",
            "--model_type=WGS",
            f"--ref={library_paths.REF_FASTA}",
            f"--reads={input_bam}",
            f"--output_vcf={output_vcf}",
            f"--output_gvcf={output_gvcf}",
            f"--num_shards={AppConfig.MAX_THREADS_PER_JOB}",
            f"--intermediate_results_dir={caller_config[VariantCallingKeys.OUTPUT_DIR]}/tmpdir/",
        ]

        bed_file = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.BED_FILE])
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )
        if bed_file is not None:
            command.append(f"--regions={bed_file}")

        return command

    @classmethod
    def create_parabricks_deepvariant_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        output_vcf = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT],
        )
        output_gvcf = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            caller_config[VariantCallingKeys.GVCF_OUTPUT],
        )

        command = [
            "pbrun",
            "deepvariant",
            f"--ref={library_paths.REF_FASTA}",
            f"--in-bam={input_bam}",
            f"--out-variants={output_vcf}",
        ]

        bed_file = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.BED_FILE])
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )
        if bed_file is not None:
            command.append(f"--interval-file={bed_file}")

        return command

    @classmethod
    def call_variants(cls, caller_config: Dict, device: str = "cpu"):
        library_paths = LibraryPaths()
        os.makedirs(caller_config[VariantCallingKeys.OUTPUT_DIR], exist_ok=True)

        workdir = os.path.abspath(caller_config[VariantCallingKeys.WORKDIR])
        input_dir = convert_to_absolute_path(
            os.path.dirname(caller_config[VariantCallingKeys.GERMLINE_INPUT])
        )

        docker_runner = DockerRunner(device=device)

        if device == "gpu":
            deepvariant_command = cls.create_parabricks_deepvariant_command(
                caller_config, library_paths
            )
            docker_runner.run(
                DockerImages.PARABRICKS,
                " ".join(deepvariant_command),
                workdir=str(workdir),
                paths_to_bind=[
                    {"source": input_dir, "target": input_dir},
                ],
            )
        else:
            deepvariant_command = cls.create_run_deepvariant_command(
                caller_config, library_paths
            )
            docker_runner.run(
                DockerImages.DEEPVARIANT,
                " ".join(deepvariant_command),
                workdir=str(workdir),
                paths_to_bind=[
                    {"source": input_dir, "target": input_dir},
                ],
            )
