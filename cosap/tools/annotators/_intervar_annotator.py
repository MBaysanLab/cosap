from subprocess import run
from typing import Dict, List

from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ..._utils import join_paths
from ._annotators import _Annotatable, _Annotator


class IntervarAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:
        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]

        intervar_db = join_paths(library_paths.INTERVAR, "intervardb")
        annovar_db = join_paths(library_paths.ANNOVAR, "humandb38")
        annotate_variation = join_paths(library_paths.ANNOVAR, "annotate_variation.pl")
        table_annovar = join_paths(library_paths.ANNOVAR, "table_annovar.pl")
        convert2annovar = join_paths(library_paths.ANNOVAR, "convert2annovar.pl")

        input_file_type = annotator_config[AnnotatorKeys.INPUT_TYPE]
        if input_file_type == "vcf":
            input_file = cls.chr_filter_vcf(input_vcf)
        elif input_file_type.lower() == "avinput":
            input_file = input_vcf

        command = [
            join_paths(library_paths.INTERVAR, "Intervar.py"),
            "-b",
            "hg38",
            "-i",
            input_file,
            f"--input_type={input_file_type}",
            "-o",
            output_vcf,
            f"--database_intervar={intervar_db}",
            "-d",
            annovar_db,
            f"--annotate_variation={annotate_variation}",
            f"--table_annovar={table_annovar}",
            f"--convert2annovar={convert2annovar}",
        ]
        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()
        workdir = annotator_config[AnnotatorKeys.OUTPUT_DIR]

        intervar_command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )

        cls.create_output_dir(annotator_config, workdir=workdir)
        run(intervar_command, cwd=workdir)
