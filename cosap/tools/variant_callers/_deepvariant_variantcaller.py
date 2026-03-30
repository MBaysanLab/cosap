import os
from typing import Dict

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ..._utils import convert_to_absolute_path
from ...runners.docker_runner import DockerRunner
from ._variantcallers import _Callable, _VariantCaller


class DeepVariantVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def create_run_deepvariant_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        # Absolute path — Docker workdir is OUTPUT_DIR, so relative BAM paths
        # would not resolve correctly.
        input_bam = convert_to_absolute_path(
            caller_config[VariantCallingKeys.GERMLINE_INPUT]
        )
        # Bare filenames — resolved relative to Docker workdir (OUTPUT_DIR).
        output_vcf = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        gvcf_mode = caller_config.get(VariantCallingKeys.OUTPUT_TYPE, "VCF") == "GVCF"
        tmpdir = os.path.abspath(
            os.path.join(caller_config[VariantCallingKeys.OUTPUT_DIR], "tmpdir")
        )

        command = [
            "run_deepvariant",
            f"--model_type={caller_config[VariantCallingKeys.MODEL_TYPE]}",
            f"--ref={library_paths.REF_FASTA}",
            f"--reads={input_bam}",
            f"--output_vcf={output_vcf}",
            f"--num_shards={AppConfig.MAX_THREADS_PER_JOB}",
            f"--intermediate_results_dir={tmpdir}/",
        ]

        if gvcf_mode:
            output_gvcf = caller_config[VariantCallingKeys.GVCF_OUTPUT]
            command.append(f"--output_gvcf={output_gvcf}")

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
        # Absolute path — Docker workdir is OUTPUT_DIR, so relative BAM paths
        # would not resolve correctly.
        input_bam = convert_to_absolute_path(
            caller_config[VariantCallingKeys.GERMLINE_INPUT]
        )
        # Bare filenames — resolved relative to Docker workdir (OUTPUT_DIR).
        output_vcf = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        gvcf_mode = caller_config.get(VariantCallingKeys.OUTPUT_TYPE, "VCF") == "GVCF"

        command = [
            "pbrun",
            "deepvariant",
            f"--ref={library_paths.REF_FASTA}",
            f"--in-bam={input_bam}",
            f"--out-variants={output_vcf}",
        ]

        if gvcf_mode:
            output_gvcf = caller_config[VariantCallingKeys.GVCF_OUTPUT]
            command.append(f"--gvcf")
            command.append(f"--out-variants={output_gvcf}")

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

        workdir = os.path.abspath(caller_config[VariantCallingKeys.OUTPUT_DIR])
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
                deepvariant_command,
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
                deepvariant_command,
                workdir=str(workdir),
                paths_to_bind=[
                    {"source": input_dir, "target": input_dir},
                ],
            )
