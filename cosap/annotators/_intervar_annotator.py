from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import AnnotatorKeys
from .._utils import join_paths
from ._annotators import _Annotatable, _Annotator


class IntervarAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:

        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        annotate_variation = join_paths(library_paths.ANNOVAR, "humandb38")
        annovar_db = join_paths(library_paths.ANNOVAR, "annotate_variation.pl")

        command = [
            library_paths.INTERVAR,
            "-b",
            "hg38",
            "-i",
            input_vcf,
            "--input_type=VCF",
            "-o",
            output_vcf,
            join_paths(library_paths.INTERVAR, "intervardb"),
            f"--annotate_variation={annotate_variation}",
            "-d",
            annovar_db,
            "--skip_annovar",
        ]

        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()

        intervar_command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )
        run(intervar_command)
