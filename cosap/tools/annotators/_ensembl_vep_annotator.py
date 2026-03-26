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
                "use_file_prefix": True,
            },
            {
                "name": "REVEL",
                "file": "REVEL_grch38.tsv.gz",
                "params": None,
                "use_file_prefix": True,
            },
            {
                "name": "SpliceAI",
                "files": {
                    "snv": "spliceai_scores.raw.snv.hg38.vcf.gz",
                    "indel": "spliceai_scores.raw.indel.hg38.vcf.gz",
                },
                "params": None,
                "use_file_prefix": True,
            },
            {
                "name": "CADD",
                "files": {
                    "snv": "CADD_whole_genome_SNVs.tsv.gz",
                    "indels": "CADD_gnomad.genomes.r4.0.indel.tsv.gz",
                },
                "params": None,
                "use_file_prefix": True,
            },
            {
                "name": "BayesDel",
                "file": "BayesDel_hg38_all_scores.txt.gz",
                "params": None,
                "use_file_prefix": True,
            },
            {
                "name": "PrimateAI",
                "file": "PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz",
                "params": None,
                "use_file_prefix": False,  # PrimateAI doesn't use file= prefix
            },
            {
                "name": "InterVar",
                "file": "intervar_db",
                "params": None,
                "use_file_prefix": True,
            },
        ]

        for plugin in plugins:
            use_file_prefix = plugin.get("use_file_prefix", True)

            # Handle single file plugins
            if "file" in plugin:
                plugin_file = join_paths(
                    library_paths.ENSEMBL_VEP, "plugin_data", plugin["file"]
                )
                if os.path.exists(plugin_file):
                    # Build plugin configuration based on whether to use file= prefix
                    if use_file_prefix:
                        plugin_config = f"{plugin['name']},file={plugin_file}"
                    else:
                        plugin_config = f"{plugin['name']},{plugin_file}"

                    if plugin["params"]:
                        plugin_config += f",{plugin['params']}"
                    command.extend(["--plugin", plugin_config])
                else:
                    print(f"Plugin file not found for {plugin['name']}: {plugin_file}")

            # Handle multiple files plugins (like SpliceAI, CADD)
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
                    if use_file_prefix:
                        file_params = ",".join(
                            [f"{key}={path}" for key, path in file_paths.items()]
                        )
                    else:
                        # If no file prefix needed, just join the paths
                        file_params = ",".join(file_paths.values())

                    plugin_config = f"{plugin['name']},{file_params}"

                    if plugin["params"]:
                        plugin_config += f",{plugin['params']}"

                    command.extend(["--plugin", plugin_config])
                else:
                    missing_files = [
                        name
                        for key, name in plugin["files"].items()
                        if key not in file_paths
                    ]
                    print(
                        f"Plugin files not found for {plugin['name']}: {', '.join(missing_files)}"
                    )

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
            "--dir_plugins",
            join_paths(library_paths.ENSEMBL_VEP, "Plugins"),
            "--cache",
            "--refseq",
            "--offline",
            "-e",
            "--fork",
            str(app_config.MAX_THREADS_PER_JOB),
            "--force_overwrite",
            "--show_ref_allele",
            "--flag_pick",
            "--cache_version",
            app_config.VEP_CACHE_VERSION,
            "--no_stats",
            "--vcf",
            "--buffer_size",
            "10000",
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
        run(["chmod", "-R", "a+rwx", output_dir], check=True)

        # Change the read/write permissions of the input file
        run(["chmod", "a+rwx", annotator_config[AnnotatorKeys.INPUT]], check=True)

        runner = DockerRunner()
        runner.run(
            DockerImages.ENSEMBL_VEP,
            cls.create_command(library_paths, app_config, annotator_config),
            workdir=output_dir,
        )
