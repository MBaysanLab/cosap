import os
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ..._utils import convert_to_absolute_path, join_paths
from ._variantcallers import _Callable, _VariantCaller


class VarNetVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_filter_command(
        cls, caller_config: Dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> List:
        varnet_filter = join_paths(library_paths.VARNET, "filter.py")
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

        sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.TUMOR_SAMPLE_NAME
        ]
        output_file = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        output_dir, _ = os.path.split(output_file)

        bed_file = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.BED_FILE])
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )

        command = [
            "python",
            varnet_filter,
            "--sample_name",
            sample_name,
            "--normal_bam",
            germline_bam,
            "--tumor_bam",
            tumor_bam,
            "--processes",
            "2",
            "--output_dir",
            output_dir,
            "--reference",
            library_paths.REF_FASTA,
            "-snv",
        ]
        if bed_file is not None:
            command.extend(["--region_bed", bed_file])
        return command

    @classmethod
    def _create_predict_command(
        cls, caller_config: Dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> List:
        varnet_predict = join_paths(library_paths.VARNET, "predict.py")

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

        sample_name = caller_config[VariantCallingKeys.PARAMS][
            VariantCallingKeys.TUMOR_SAMPLE_NAME
        ]

        output_file = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_dir, _ = os.path.split(output_file)
        command = [
            "python",
            varnet_predict,
            "--sample_name",
            sample_name,
            "--normal_bam",
            germline_bam,
            "--tumor_bam",
            tumor_bam,
            "--processes",
            "2",
            "--output_dir",
            output_dir,
            "--output_vcf",
            output_file,
            "--reference",
            library_paths.REF_FASTA,
            "-snv",
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict, device: str = "cpu"):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        varnet_filter_command = cls._create_filter_command(
            caller_config=caller_config,
            library_paths=library_paths,
            app_config=app_config,
        )
        varnet_predict_command = cls._create_predict_command(
            caller_config=caller_config,
            library_paths=library_paths,
            app_config=app_config,
        )

        workdir = caller_config[VariantCallingKeys.OUTPUT_DIR]
        run(varnet_filter_command, cwd=workdir)
        run(varnet_predict_command, cwd=workdir)
