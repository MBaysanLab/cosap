import os
from pathlib import Path
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ..._utils import join_paths
from ...pipeline_runner.runners import DockerRunner
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
            "--fa",
            library_paths.REF_FASTA,
            "-i",
            input_vcf,
            "-o",
            output_vcf,
            "--dir",
            library_paths.ENSEMBL_VEP,
            "--dir_cache",
            join_paths(library_paths.ENSEMBL_VEP, "cache"),
            "--fasta",
            library_paths.REF_FASTA,
            "--cache",
            "--offline",
            "--tab",
            "-e",
            "--fork",
            str(app_config.MAX_THREADS_PER_JOB),
            "--force_overwrite",
            "--show_ref_allele",
        ]
        return command

    @classmethod
    def annotate(cls, annotator_config: Dict, workdir: str = None):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        output_dir = cls.create_output_dir(annotator_config, workdir=workdir)

        runner = DockerRunner()
        runner.run(
            DockerImages.ENSEMBL_VEP,
            " ".join(cls.create_command(library_paths, app_config, annotator_config)),
            workdir=str(Path(output_dir).parent.parent) if not workdir else workdir,
        )
