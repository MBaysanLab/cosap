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

        annovar_db = join_paths(library_paths.ANNOVAR, "humandb38")
        intervar_db = join_paths(library_paths.INTERVAR, "intervardb")
        table_annovar = annotate_variation = annovar_db = join_paths(
            library_paths.ANNOVAR, "table_annovar.pl"
        )

        command = [
            join_paths(library_paths.INTERVAR, "Intervar.py"),
            "-b",
            "hg38",
            "-i",
            input_vcf,
            "--input_type=VCF",
            "-o",
            output_vcf,
            f"--database_intervar={intervar_db}",
            "-d",
            annovar_db,
            "--table-annovar",
            table_annovar
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
