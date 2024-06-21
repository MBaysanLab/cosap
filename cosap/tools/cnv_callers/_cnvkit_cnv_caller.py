import glob
import os
from subprocess import run

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import CNVCallingKeys
from ._cnv_callers import _CNVCaller


class CNVKit(_CNVCaller):
    @classmethod
    def _create_cnvkit_batch_command(
        cls, caller_config: dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> list:
        normal_input = caller_config[CNVCallingKeys.NORMAL_INPUT]
        tumor_input = caller_config[CNVCallingKeys.TUMOR_INPUT]
        bed_file = (
            caller_config[CNVCallingKeys.BED_FILE]
            if CNVCallingKeys.BED_FILE in caller_config
            else None
        )
        output_dir = caller_config[CNVCallingKeys.OUTPUT_DIR]

        command = [
            "cnvkit.py",
            "batch",
            *tumor_input,
            "--normal",
            *normal_input,
            "--targets" if bed_file is None else f"--targets {bed_file}",
            "--annotate",
            library_paths.CNVKIT_ANNOTATION,
            "--fasta",
            library_paths.REF_FASTA,
            "--access",
            library_paths.CNVKIT_ACCESS,
            "--output-reference",
            library_paths.CNVKIT_REFERENCE,
            "--output-dir",
            output_dir,
        ]
        return command

    @classmethod
    def _create_cnvkit_genemetrics_command(
        cls, caller_config: dict, library_paths: LibraryPaths, app_config: AppConfig
    ):
        input_cnr = glob.glob(
            os.path.join(caller_config[CNVCallingKeys.OUTPUT_DIR], "*.cnr")
        )
        output = caller_config[CNVCallingKeys.OUTPUT]

        command = ["cnvkit.py", "genemetrics", *input_cnr, ">", output]

        return command

    @classmethod
    def call(cls, caller_config: dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        batch_command = cls._create_cnvkit_batch_command(
            caller_config, library_paths, app_config
        )
        run(batch_command)
        genemetrics_command = cls._create_cnvkit_genemetrics_command(
            caller_config, library_paths, app_config
        )
        run(genemetrics_command)
