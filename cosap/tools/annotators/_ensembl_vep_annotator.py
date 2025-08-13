from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import AnnotatorKeys
from ..._utils import join_paths
from ...runners.docker_runner import DockerRunner
from ._annotators import _Annotatable, _Annotator


class VepAnnotator(_Annotatable, _Annotator):
    @classmethod
    def create_command(
        cls, library_paths: LibraryPaths, app_config: AppConfig, annotator_config: Dict
    ) -> List:
        input_vcf = annotator_config[AnnotatorKeys.INPUT]
        output_vcf = annotator_config[AnnotatorKeys.OUTPUT]
        input_format = annotator_config[AnnotatorKeys.INPUT_TYPE]

        command = [
            "vep",
            "-v",
            "--fa",
            library_paths.REF_FASTA,
            "-i",
            input_vcf,
            "--format",
            input_format,
            "-o",
            output_vcf,
            "--dir",
            library_paths.ENSEMBL_VEP,
            "--dir_cache",
            join_paths(library_paths.ENSEMBL_VEP, "cache"),
            "--cache",
            "--offline",
            "-e",
            "--fork",
            str(app_config.MAX_THREADS_PER_JOB),
            "--force_overwrite",
            "--show_ref_allele",
            "--flag_pick",
            "--hgvsg",
            "--hgvsg_use_accession",
            "--cache_version",
            app_config.VEP_CACHE_VERSION,
            "--plugin",
            f"AlphaMissense,file={join_paths(library_paths.ENSEMBL_VEP, 'AlphaMissense_hg38.tsv.gz')},pathogenicity_cutoff=0.5",
            "--no_stats",
            "--vcf",
        ]
        return command

    @classmethod
    def annotate(cls, annotator_config: Dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()

        output_dir = cls.create_output_dir(annotator_config)

        # Change the read/write permissions of the output directory
        run(["chmod", "-R", "a+rwx", output_dir])

        # Change the read/write permissions of the input file
        run(["chmod", "a+rwx", annotator_config[AnnotatorKeys.INPUT]])

        runner = DockerRunner()
        runner.run(
            DockerImages.ENSEMBL_VEP,
            " ".join(cls.create_command(library_paths, app_config, annotator_config)),
            workdir=output_dir,
        )
