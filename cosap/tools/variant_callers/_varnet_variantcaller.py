import os
import shutil
from pathlib import Path

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ..._utils import current_available_system_memory
from ...pipeline_runner.runners import DockerRunner
from ._variantcallers import _Callable, _VariantCaller


class VarNetVariantCaller(_Callable, _VariantCaller):
    # https://github.com/skandlab/VarNet#requirements
    FILTER_THREAD_MEMORY_GB = 6
    PREDICT_THREAD_MEMORY_GB = 10

    GIGABYTE = 1_000_000_000
    FILTER_THREAD_MEMORY = FILTER_THREAD_MEMORY_GB * GIGABYTE
    PREDICT_THREAD_MEMORY = PREDICT_THREAD_MEMORY_GB * GIGABYTE

    valid_variant_types = ["snv", "indel"]

    @classmethod
    def _get_tumor_sample_name(cls, caller_config: dict):
        tumor_sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.TUMOR_SAMPLE_NAME
        ]
        return tumor_sample_name

    @classmethod
    def _get_thread_count(cls, app_config: AppConfig, memory_per_thread: int) -> int:
        available_memory = 0.9 * current_available_system_memory()  # Tolerate minor fluctuations

        maximum_acceptable_thread_count = available_memory // memory_per_thread

        target_thread_count = min(maximum_acceptable_thread_count, app_config.MAX_THREADS_PER_JOB)

        if target_thread_count == 0:
            available_memory_GB = int(available_memory // cls.GIGABYTE)
            memory_per_thread_GB = int(memory_per_thread // cls.GIGABYTE)
            raise ResourceWarning(f"Insufficient system memory. "
                                  f"Minimum Required: {memory_per_thread_GB}GB, "
                                  f"Current Available: {available_memory_GB}GB")

        return int(target_thread_count)
    
    @classmethod
    def _populate_common_arguments(
        cls, caller_config: dict, library_paths: LibraryPaths, app_config: AppConfig, variant_type: str
    ) -> list:
        germline_bam = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.GERMLINE_INPUT])
            if VariantCallingKeys.GERMLINE_INPUT in caller_config.keys()
            else None
        )
        tumor_bam = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.TUMOR_INPUT])
            if VariantCallingKeys.TUMOR_INPUT in caller_config.keys()
            else None
        )

        tumor_sample_name = cls._get_tumor_sample_name(caller_config)

        output_file = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_dir, _ = os.path.split(output_file)

        if variant_type not in cls.valid_variant_types:
            raise ValueError(f"Unexpected variant_type: Expected one of {cls.valid_variant_types}, got {variant_type}.")

        common_arguments = [
            "--sample_name", tumor_sample_name,
            "--normal_bam", germline_bam,
            "--tumor_bam", tumor_bam,
            "--output_dir", output_dir,
            "--reference", library_paths.REF_FASTA,
            f"-{variant_type}"
        ]
        
        return common_arguments
    
    @classmethod
    def _create_filter_command(
        cls, caller_config: dict, library_paths: LibraryPaths, app_config: AppConfig, variant_type: str = "snv"
    ) -> list:
        varnet_filter_script = "/VarNet/filter.py"  # Location within Docker container

        bed_file = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.BED_FILE])
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )

        thread_count = cls._get_thread_count(app_config, memory_per_thread=cls.FILTER_THREAD_MEMORY)

        common_arguments = cls._populate_common_arguments(caller_config, library_paths, app_config, variant_type)

        command = [
            "python", varnet_filter_script,
            *common_arguments,
            "--processes", str(thread_count)
        ]
        if bed_file is not None:
            command.extend(["--region_bed", bed_file])
        return command

    @classmethod
    def _create_predict_command(
        cls, caller_config: dict, library_paths: LibraryPaths, app_config: AppConfig, variant_type: str = "snv"
    ) -> list:
        varnet_predict_script = "/VarNet/predict.py"  # Location within Docker container

        thread_count = cls._get_thread_count(app_config, memory_per_thread=cls.PREDICT_THREAD_MEMORY)

        common_arguments = cls._populate_common_arguments(caller_config, library_paths, app_config, variant_type)

        command = [
            "python", varnet_predict_script,
            *common_arguments,
            "--processes", str(thread_count)
        ]
        return command
    
    @classmethod
    def _create_rename_output_vcf_command(cls, caller_config: dict, variant_type: str = "snv"):
            
            tumor_sample_name = cls._get_tumor_sample_name(caller_config)

            varnet_output_vcf = f"VCF/varnet/{tumor_sample_name}/{tumor_sample_name}.vcf"
            varnet_renamed_output_vcf =  f"VCF/varnet/{tumor_sample_name}/{tumor_sample_name}.{variant_type}.vcf"

            return ["mv", varnet_output_vcf, varnet_renamed_output_vcf]

    @classmethod
    def _move_vcf(
        cls, caller_config: dict, variant_type: str = "snv"
    ) -> list:
        tumor_sample_name = cls._get_tumor_sample_name(caller_config)

        source_vcf_path = f"VCF/varnet/{tumor_sample_name}/{tumor_sample_name}.{variant_type}.vcf"

        if variant_type not in cls.valid_variant_types:
            raise ValueError(f"Unexpected variant_type: Expected one of {cls.valid_variant_types}, got {variant_type}.")

        if variant_type == "snv":
            destination_vcf_path = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        
        if variant_type == "indel":
            destination_vcf_path = caller_config[VariantCallingKeys.INDEL_OUTPUT]
            
        with open(source_vcf_path, "rb") as source_vcf:
            with open(destination_vcf_path, "wb") as destination_vcf:
                shutil.copyfileobj(source_vcf, destination_vcf)

    @classmethod
    def call_variants(cls, caller_config: dict, device: str = "cpu"):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        output_dir = os.path.abspath(
            os.path.dirname(caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT])
        )
        workdir = Path(output_dir).parent.parent

        for variant_type in cls.valid_variant_types:
            varnet_filter_command = cls._create_filter_command(
                caller_config=caller_config,
                library_paths=library_paths,
                app_config=app_config,
                variant_type=variant_type
            )
            varnet_predict_command = cls._create_predict_command(
                caller_config=caller_config,
                library_paths=library_paths,
                app_config=app_config,
                variant_type=variant_type
            )

            varnet_output_rename_command = cls._create_rename_output_vcf_command(caller_config, variant_type)

            docker_runner = DockerRunner(device=device)
            docker_runner.run(
                DockerImages.VARNET,
                " ".join(varnet_filter_command),
                workdir=str(workdir)
            )

            docker_runner.run(
                DockerImages.VARNET,
                " ".join(varnet_predict_command),
                workdir=str(workdir)
            )

            docker_runner.run(
                DockerImages.VARNET,
                " ".join(varnet_output_rename_command),
                workdir=str(workdir)
            )

            cls._move_vcf(
                caller_config=caller_config,
                variant_type=variant_type
            )
