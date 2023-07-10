from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import AnnotatorKeys
from .._utils import join_paths
from ._annotators import _Annotatable, _Annotator
import os


class AnnotSVAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, annotator_config: Dict
    ) -> List:

        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_file = annotator_config[AnnotatorKeys.OUTPUT]


        command = [
            library_paths.ANNOTSV,
            "-SVinputFile",
            input_vcf,
            "-outputFile",
            output_file,
            "-outputDir",
            os.getcwd()
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
