from subprocess import run
from typing import Dict, List

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import AnnotatorKeys
from ._annotators import _Annotatable, _Annotator
from .._utils import join_paths


class AnnovarAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_av_input_command(cls,library_paths: LibraryPaths, annotator_config: Dict):
        command = [
            join_paths(library_paths.ANNOVAR, "convert2annovar.pl"),
            "-format",
            "vcf4",
            annotator_config[AnnotatorKeys.OUTPUT],
            "-outfile",
            annotator_config[AnnotatorKeys.AVOUTPUT]
        ]
        return command

    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:

        input_vcf = annotator_config[AnnotatorKeys.AVOUTPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        command = [
            join_paths(library_paths.ANNOVAR, "annotate_variation.pl"),
            "-geneanno",
            "-dbtype",
            "refGene",
            "-buildver",
            "hg38",
            input_vcf,
            join_paths(library_paths.ANNOVAR, "humandb38"),
            "--outfile",
            output_vcf
        ]

        return command


    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()

        create_input_command = cls.create_av_input_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )
        annovar_command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )
        run(create_input_command)
        run(annovar_command)
