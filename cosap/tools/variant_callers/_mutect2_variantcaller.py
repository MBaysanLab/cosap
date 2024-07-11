import os
from itertools import chain, repeat
from pathlib import Path
from subprocess import run

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ...memory_handler import MemoryHandler
from ...pipeline_runner.runners import DockerRunner
from ...scatter_gather import ScatterGather
from ._variantcallers import _Callable, _VariantCaller


class Mutect2VariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls,
        caller_config: dict,
        library_paths: LibraryPaths,
        memory_handler: MemoryHandler,
    ) -> list:
        MAX_MEMORY_IN_GB = int(AppConfig.MAX_MEMORY_PER_JOBS // (1024**3))

        germline_bam = None
        if VariantCallingKeys.GERMLINE_INPUT in caller_config.keys():
            germline_input = caller_config[VariantCallingKeys.GERMLINE_INPUT]
            germline_bam = memory_handler.get_bam_path(germline_input)

        tumor_bam = None
        if VariantCallingKeys.TUMOR_INPUT in caller_config.keys():
            tumor_input = caller_config[VariantCallingKeys.TUMOR_INPUT]
            tumor_bam = memory_handler.get_bam_path(tumor_input)

        output_name = memory_handler.get_output_path(
            caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        )

        ref_fasta = library_paths.REF_FASTA

        command = [
            "gatk",
            "--java-options",
            f"-Xmx{MAX_MEMORY_IN_GB}G -XX:+UseParallelGC -XX:ParallelGCThreads=1",
            "Mutect2",
            "-R",
            ref_fasta,
            "-I",
            tumor_bam,
            "-O",
            output_name,
            "--native-pair-hmm-threads",
            "2",
        ]

        if germline_bam is not None:
            germline_sample_name = caller_config[VariantCallingKeys.PARAMS][
                VariantCallingKeys.GERMLINE_SAMPLE_NAME
            ]
            command.extend(
                [
                    "-I",
                    germline_bam,
                    "-normal",
                    germline_sample_name,
                ]
            )

        bed_file = (
            caller_config[VariantCallingKeys.BED_FILE]
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )
        if bed_file is not None:
            command.extend(["--intervals", bed_file])

        return command

    @classmethod
    def _create_get_snp_variants_command(
        cls, caller_config: dict, library_paths: LibraryPaths
    ) -> list:
        input_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.SNP_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
            "-V",
            input_name,
            "--exclude-filtered",
            "--select-type-to-include",
            "SNP",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_indel_variants_command(
        cls, caller_config: dict, library_paths: LibraryPaths
    ) -> list:
        input_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.INDEL_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
            "-V",
            input_name,
            "--exclude-filtered",
            "--select-type-to-include",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_other_variants_command(
        cls, caller_config: dict, library_paths: LibraryPaths
    ) -> list:
        input_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.OTHER_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
            "-V",
            input_name,
            "--exclude-filtered",
            "--select-type-to-exclude",
            "SNP",
            "--select-type-to-exclude",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _gather_mutect_stats_command(cls, configs: list, output) -> list:
        stats_files = []
        for cfg in configs:
            stats_file = f"{cfg[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]}.stats"
            if os.path.exists(stats_file):
                stats_files.append(stats_file)

        command = ["gatk", "MergeMutectStats", "-O", output]
        command.extend(list(chain(*zip(repeat("--stats"), stats_files))))
        return command

    @classmethod
    def _create_parabricks_mutectcaller_command(
        cls, caller_config: dict, library_paths: LibraryPaths
    ) -> list:

        germline_bam = None
        if VariantCallingKeys.GERMLINE_INPUT in caller_config.keys():
            germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        tumor_bam = None
        if VariantCallingKeys.TUMOR_INPUT in caller_config.keys():
            tumor_bam = caller_config[VariantCallingKeys.TUMOR_INPUT]

        output_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]

        command = [
            "pbrun",
            "mutectcaller",
            "--ref",
            library_paths.REF_FASTA,
            "--tumor-name",
            caller_config[VariantCallingKeys.PARAMS][
                VariantCallingKeys.TUMOR_SAMPLE_NAME
            ],
            "--in-tumor-bam",
            tumor_bam,
            "--out-vcf",
            output_name,
        ]
        if germline_bam is not None:
            command.extend(
                [
                    "--normal-name",
                    caller_config[VariantCallingKeys.PARAMS][
                        VariantCallingKeys.GERMLINE_SAMPLE_NAME
                    ],
                    "--in-normal-bam",
                    germline_bam,
                ]
            )

        return command

    @classmethod
    def _filter_mutect_calls(
        cls, caller_config: dict, library_paths: LibraryPaths
    ) -> list:
        input_name = caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "FilterMutectCalls",
            "-V",
            input_name,
            "-R",
            library_paths.REF_FASTA,
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def call_variants(cls, caller_config: dict, device: str = "cpu"):
        library_paths = LibraryPaths()

        if device == "cpu":
            bed_file = (
                caller_config[VariantCallingKeys.BED_FILE]
                if VariantCallingKeys.BED_FILE in caller_config.keys()
                else None
            )
            splitted_configs = ScatterGather.split_variantcaller_configs(
                caller_config, bed_file=bed_file
            )
            print(f"Splitting {len(splitted_configs)} configs")
            with MemoryHandler() as memory_handler:
                scattered_commands = [
                    cls._create_run_command(
                        caller_config=cfg,
                        library_paths=library_paths,
                        memory_handler=memory_handler,
                    )
                    for cfg in splitted_configs
                ]

                ScatterGather.run_parallel(run, scattered_commands)

                stats_output = f"{caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]}.stats"
                gather_stats_command = cls._gather_mutect_stats_command(
                    splitted_configs, stats_output
                )
            ScatterGather.gather_vcfs(
                splitted_configs,
                output_path=caller_config[
                    VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT
                ],
            )

            run(gather_stats_command)

            ScatterGather.clean_temp_files(caller_config[VariantCallingKeys.OUTPUT_DIR])

        elif device == "gpu":

            print("Running Mutect2 on GPU")
            command = cls._create_parabricks_mutectcaller_command(
                caller_config=caller_config, library_paths=library_paths
            )
            output_dir = os.path.abspath(
                os.path.dirname(
                    caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
                )
            )
            os.makedirs(output_dir, exist_ok=True)

            runner = DockerRunner(device=device)
            runner.run(
                image=DockerImages.PARABRICKS,
                command=" ".join(command),
                workdir=str(Path(output_dir).parent.parent),
            )

        filter_command = cls._filter_mutect_calls(
            caller_config=caller_config, library_paths=library_paths
        )
        get_snp_command = cls._create_get_snp_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )
        get_indel_command = cls._create_get_indel_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )
        get_other_variants_command = cls._create_get_other_variants_command(
            caller_config=caller_config, library_paths=library_paths
        )

        workdir = caller_config[VariantCallingKeys.OUTPUT_DIR]
        run(filter_command, cwd=workdir)
        run(get_snp_command, cwd=workdir)
        run(get_indel_command, cwd=workdir)
        run(get_other_variants_command, cwd=workdir)
