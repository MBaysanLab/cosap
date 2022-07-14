from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import AnnotatorKeys
from ._annotators import _Annotatable, _Annotator


class VepAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:

        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        command = [
            library_paths.ENSEMBL_VEP,
            "-i",
            input_vcf,
            "-o",
            output_vcf,
            "--cache",
            "--offline",
            "--vcf",
            "-e"
        ]
        print(" ".join(command))
        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()

        command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )
        run(command)
