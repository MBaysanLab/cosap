import os
from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import VariantCallingKeys
from ._variantcallers import _Callable, _VariantCaller
from .._containers import DockerContainers
import docker
from .._utils import check_if_running_in_docker
from pathlib import Path

class DeepVariantVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def create_run_deepvariant_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ):
        
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]
        output_vcf = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_gvcf = caller_config[VariantCallingKeys.GVCF_OUTPUT]

        command = [
            "run_deepvariant",
            "--model_type=WGS",
            f"--ref={library_paths.REF_FASTA}",
            f"--reads={input_bam}",
            f"--output_vcf={output_vcf}",
            f"--output_gvcf={output_gvcf}",
            f"--num_shards={AppConfig.MAX_THREADS_PER_JOB}",
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()
        docker_client = docker.from_env()
        os.makedirs(caller_config[VariantCallingKeys.OUTPUT_DIR], exist_ok=True)

        # If running in docker mount volumes from host to container
        if check_if_running_in_docker():
            hostname = os.getenv("HOSTNAME")
            docker_client.containers.run(
                image=DockerContainers.DEEPVARIANT,
                working_dir=docker_client.containers.get(hostname).attrs["Config"]["WorkingDir"],
                command=" ".join(cls.create_run_deepvariant_command(caller_config, library_paths)),
                volumes_from=[hostname],
                remove=True,
                detach=False,
            )
        else:
            # If not running in docker, mount the required paths to the container
            input_dir = os.path.abspath(os.path.dirname(caller_config[VariantCallingKeys.GERMLINE_INPUT]))
            output_dir = os.path.abspath(os.path.dirname(caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]))
            library_path = AppConfig.LIBRARY_PATH
            container = docker_client.containers.run(
                image=DockerContainers.DEEPVARIANT,
                command=" ".join(cls.create_run_deepvariant_command(caller_config, library_paths)),
                working_dir=str(Path(output_dir).parent.parent),
                volumes={
                    library_path: {"bind": library_path, "mode": "ro"},
                    input_dir: {"bind": input_dir, "mode": "ro"},
                    output_dir: {"bind": output_dir, "mode": "rw"},
                },
                remove=True,
                detach=False,
                restart_policy={"Name": "no"},
            )


