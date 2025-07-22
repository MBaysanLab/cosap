import os
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ..._utils import join_paths
from ._annotators import _Annotatable, _Annotator


class AnnovarAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_av_input_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ):
        avinput_filename = annotator_config[AnnotatorKeys.AVOUTPUT].split(".")[0]
        command = [
            join_paths(library_paths.ANNOVAR, "convert2annovar.pl"),
            "-format",
            "vcf4",
            "-allsample",
            annotator_config[AnnotatorKeys.INPUT],
            "-outfile",
            avinput_filename,
        ]
        return command

    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:
        avinput = annotator_config[AnnotatorKeys.AVOUTPUT]
        output_filename = annotator_config[AnnotatorKeys.OUTPUT].split(".")[0]

        command = [
            join_paths(library_paths.ANNOVAR, "table_annovar.pl"),
            avinput,
            join_paths(library_paths.ANNOVAR, "humandb38"),
            "--buildver",
            "hg38",
            "--out",
            output_filename,
            "--remove",
            "--protocol",
            "ensGene,gnomad40,icgc28,avsnp150,dbnsfp42a,"
            "clinvar_20221231,intervar_20180118",
            "--operation",
            "g,f,f,f,f,f,f",
            "--nastring",
            ".",
        ]
        return command

    @classmethod
    def _rename_annovar_output(self, annotator_config: Dict):
        annotation_output = annotator_config[AnnotatorKeys.OUTPUT]
        output_filename = annotator_config[AnnotatorKeys.OUTPUT].split(".")[0]
        os.rename(f"{output_filename}.hg38_multianno.txt", annotation_output)

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
        run(create_input_command, cwd=annotator_config[AnnotatorKeys.OUTPUT_DIR])
        run(annovar_command, cwd=annotator_config[AnnotatorKeys.OUTPUT_DIR])
        cls._rename_annovar_output(annotator_config=annotator_config)
