import os
from typing import Dict, List

from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ..._utils import join_paths
from ...runners.docker_runner import DockerRunner
from ._annotators import _Annotatable, _Annotator


class MetaRnnAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, annotator_config: Dict
    ) -> List:
        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        annotate_variation = join_paths(library_paths.ANNOVAR, "annotate_variation.pl")
        annovar_db = join_paths(library_paths.ANNOVAR, "humandb38")
        snv_anno = join_paths(library_paths.METARNN, "SNV_Anno")
        transcript_anno = join_paths(library_paths.METARNN, "Transcript_Anno")

        output_file = annotator_config.get(AnnotatorKeys.OUTPUT)

        command = [
            "metarnn",
            "--genome_build",
            "hg38",
            "--input_vcf",
            input_vcf,
            "--annovar_script",
            annotate_variation,
            "--annovar_db",
            annovar_db,
            "--snv_anno",
            snv_anno,
            "--transcript_anno",
            transcript_anno,
            "--output_file",
            output_file,
        ]

        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()

        output_dir = cls.create_output_dir(annotator_config)
        input_file = annotator_config[AnnotatorKeys.INPUT]
        input_dir = os.path.abspath(os.path.dirname(input_file))

        metarnn_command = cls.create_command(
            library_paths=library_paths,
            annotator_config=annotator_config,
        )

        runner = DockerRunner()

        runner.run(
            image=DockerImages.METARNN,
            command=" ".join(metarnn_command),
            workdir=output_dir,
            paths_to_bind=[
                {"source": input_dir, "target": input_dir},
            ],
        )
