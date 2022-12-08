from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import AnnotatorKeys
from .._utils import join_paths
from ._annotators import _Annotatable, _Annotator


class VepAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, annotator_config: Dict
    ) -> List:

        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        command = [
            "vep",
            "-v",
            "-i",
            input_vcf,
            "--format",
            "vcf",
            "-o",
            output_vcf,
            "--dir",
            join_paths(library_paths.ENSEMBL_VEP, "cache"),
            "--fasta",
            library_paths.REF_FASTA,
            "--cache",
            "--offline",
            "--vcf",
            "-e",
            "--fork",
            str(app_config.THREADS),
            "--plugin",
            "Phenotypes"
        ]
        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        command = cls.create_command(
            library_paths=library_paths,
            app_config=app_config,
            annotator_config=annotator_config,
        )
        run(command)
