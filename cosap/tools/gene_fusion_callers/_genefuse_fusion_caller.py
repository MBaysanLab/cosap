from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import GeneFusionCallingKeys
from ..._utils import join_paths
from ._gene_fusion_callers import _GeneFusionCaller


class GeneFuse(_GeneFusionCaller):
    @classmethod
    def _create_gene_fuse_command(
        cls, caller_config: Dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> List:
        fastq_inputs = [
            fastq for fastq in caller_config[GeneFusionCallingKeys.INPUT].values()
        ]
        output_json = caller_config[GeneFusionCallingKeys.OUTPUT]
        output_html = caller_config[GeneFusionCallingKeys.OUTPUT].replace(
            ".json", ".html"
        )

        command = [
            "genefuse",
            "-r",
            library_paths.REF_FASTA,
            "-f",
            library_paths.GENEFUSE_CANCER_GENE_LIST,
            "-1",
            fastq_inputs[0],
            "-2",
            fastq_inputs[1],
            "-j",
            output_json,
            "-h",
            output_html,
            "-t",
            str(app_config.MAX_THREADS_PER_JOB),
        ]
        return command

    @classmethod
    def call(cls, caller_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()
        command = cls._create_gene_fuse_command(
            caller_config, library_paths, app_config
        )
        run(command, cwd=caller_config[GeneFusionCallingKeys.OUTPUT_DIR])
