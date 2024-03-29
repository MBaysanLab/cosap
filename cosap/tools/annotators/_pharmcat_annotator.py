import os
from pathlib import Path
from subprocess import run
from typing import Dict, List

from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ..._utils import join_paths
from ._annotators import _Annotatable, _Annotator


class PharmcatAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_preprocess_vcf_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:
        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        command = [
            "python3",
            library_paths.PHARMCAT_PREPROCESSOR,
            "--vcf",
            input_vcf,
            "--bf",
            Path(output_vcf).stem,
        ]
        return command

    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:
        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        output_dir = os.path.abspath(
            os.path.dirname(annotator_config[AnnotatorKeys.OUTPUT])
        )

        command = [
            "java",
            "-jar",
            library_paths.PHARMCAT_JAR,
            "-vcf",
            input_vcf,
            "-bf",
            Path(output_vcf).stem,
            "-o",
            str(Path(output_dir)),
        ]

        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()

        os.makedirs(annotator_config[AnnotatorKeys.OUTPUT], exist_ok=True)

        create_preprocess_vcf_command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )
        pharmcat_command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )

        run(create_preprocess_vcf_command)
        run(pharmcat_command)
