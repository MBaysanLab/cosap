import os
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
    def _add_vep_plugins(cls, command: List[str], library_paths: LibraryPaths) -> None:
        """
        Add VEP plugins to the command if their required files exist.

        Args:
            command: The VEP command list to extend
            library_paths: Library paths configuration
        """
        # Define plugins with their file requirements and parameters
        plugins = [
            {
                "name": "AlphaMissense",
                "file": "AlphaMissense_hg38.tsv.gz",
                "params": "pathogenicity_cutoff=0.5",
            },
            {"name": "REVEL", "file": "REVEL_grch38.tsv.gz", "params": None},
            {
                "name": "SpliceAI",
                "files": {
                    "snv": "spliceai_scores.raw.snv.hg38.vcf.gz",
                    "indel": "spliceai_scores.raw.indel.hg38.vcf.gz",
                },
                "params": None,
            },
            # Add more plugins here as needed
            # {
            #     "name": "CADD",
            #     "file": "CADD_GRCh38.tsv.gz",
            #     "params": None
            # },
        ]

        for plugin in plugins:
            # Handle single file plugins
            if "file" in plugin:
                plugin_file = join_paths(
                    library_paths.ENSEMBL_VEP, "plugin_data", plugin["file"]
                )
                if os.path.exists(plugin_file):
                    plugin_config = f"{plugin['name']},file={plugin_file}"
                    if plugin["params"]:
                        plugin_config += f",{plugin['params']}"
                    command.extend(["--plugin", plugin_config])

            # Handle multiple files plugins (like SpliceAI)
            elif "files" in plugin:
                file_paths = {}
                all_files_exist = True

                # Check if all required files exist
                for file_key, file_name in plugin["files"].items():
                    file_path = join_paths(
                        library_paths.ENSEMBL_VEP, "plugin_data", file_name
                    )
                    if os.path.exists(file_path):
                        file_paths[file_key] = file_path
                    else:
                        all_files_exist = False
                        break

                if all_files_exist:
                    # Build plugin configuration with file keys
                    file_params = ",".join(
                        [f"{key}={path}" for key, path in file_paths.items()]
                    )
                    plugin_config = f"{plugin['name']},{file_params}"

                    if plugin["params"]:
                        plugin_config += f",{plugin['params']}"

                    command.extend(["--plugin", plugin_config])

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
            "--assembly",
            "GRCh38",
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
            "--no_stats",
            "--vcf",
        ]

        # Add plugins dynamically based on file availability
        cls._add_vep_plugins(command, library_paths)

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
