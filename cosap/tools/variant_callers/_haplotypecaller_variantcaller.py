import os
from pathlib import Path
from subprocess import run
from typing import Dict, List

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import VariantCallingKeys
from ..._utils import convert_to_absolute_path, join_paths
from ...memory_handler import MemoryHandler
from ...runners.runners import DockerRunner
from ...runners.runners.step_runner import run_command_parallel
from ...scatter_gather import ScatterGather
from ._variantcallers import _Callable, _VariantCaller


class HaplotypeCallerVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_run_command(
        cls,
        caller_config: Dict,
        library_paths: LibraryPaths,
        memory_handler: MemoryHandler,
    ) -> List:
        MAX_MEMORY_IN_GB = int(AppConfig.MAX_MEMORY_PER_JOB // (1024.0**3))

        germline_bam = convert_to_absolute_path(
            memory_handler.get_bam_path(
                caller_config[VariantCallingKeys.GERMLINE_INPUT]
            )
        )
        output_file = (
            caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
            if caller_config[VariantCallingKeys.OUTPUT_TYPE] == "VCF"
            else caller_config[VariantCallingKeys.GVCF_OUTPUT]
        )

        output_name = memory_handler.get_output_path(output_file)

        command = [
            "gatk",
            "--java-options",
            f"-Xmx20G",
            "HaplotypeCaller",
            "-R",
            library_paths.REF_FASTA,
            "-I",
            germline_bam,
            "-O",
            output_name,
        ]
        if caller_config[VariantCallingKeys.OUTPUT_TYPE] == "GVCF":
            command.append("--emit-ref-confidence")
            command.append("GVCF")

        bed_file = (
            convert_to_absolute_path(caller_config[VariantCallingKeys.BED_FILE])
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )
        if bed_file is not None:
            command.extend(["--intervals", bed_file])

        return command

    @classmethod
    def _create_parabricks_haplotypecaller_command(
        cls, caller_config: dict, library_paths: LibraryPaths
    ) -> list:

        germline_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        output_file = (
            caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
            if caller_config[VariantCallingKeys.OUTPUT_TYPE] == "VCF"
            else caller_config[VariantCallingKeys.GVCF_OUTPUT]
        )

        command = [
            "pbrun",
            "haplotypecaller",
            "--ref",
            library_paths.REF_FASTA,
            "--in-bam",
            germline_bam,
            "--out-variants",
            output_file,
        ]
        if caller_config[VariantCallingKeys.OUTPUT_TYPE] == "GVCF":
            command.append("--gvcf")

        return command

    @classmethod
    def _create_cnn_annotated_output_name(cls, output_name: str) -> str:
        return output_name.replace(".vcf", ".cnn_annotated.vcf")

    @classmethod
    def _create_cnnscorevariants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:
        input_name = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT],
        )
        output_name = join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            cls._create_cnn_annotated_output_name(
                caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
            ),
        )
        input_bam = caller_config[VariantCallingKeys.GERMLINE_INPUT]

        command = [
            "gatk",
            "NVScoreVariants",
            "-I",
            input_bam,
            "-V",
            input_name,
            "-R",
            library_paths.REF_FASTA,
            "-O",
            output_name,
            "-tensor-type",
            "read_tensor",
        ]

        return command
    
    @classmethod
    def _create_vcf_index_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:
        input_name = convert_to_absolute_path(join_paths(
            caller_config[VariantCallingKeys.OUTPUT_DIR],
            cls._create_cnn_annotated_output_name(
                caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
            ),
        ))

        command = [
            "gatk",
            "IndexFeatureFile",
            "-I",
            input_name,
        ]

        return command

    @classmethod
    def _create_filter_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:
        input_name = convert_to_absolute_path(
            join_paths(
                caller_config[VariantCallingKeys.OUTPUT_DIR],
                cls._create_cnn_annotated_output_name(
                    caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
                ),
            )
        )
        output_name = caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "FilterVariantTranches",
            "-V",
            input_name,
            "--resource",
            library_paths.MILLS_INDEL,
            "--resource",
            library_paths.DBSNP,
            "--resource",
            library_paths.ONE_THOUSAND_G,
            "--info-key",
            "CNN_2D",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_snp_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:
        input_name = convert_to_absolute_path(
            join_paths(
                caller_config[VariantCallingKeys.OUTPUT_DIR],
                caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT],
            )
        )
        output_name = caller_config[VariantCallingKeys.SNP_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
            "-V",
            input_name,
            "--select-type-to-include",
            "SNP",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_indel_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:
        input_name = convert_to_absolute_path(
            join_paths(
                caller_config[VariantCallingKeys.OUTPUT_DIR],
                caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT],
            )
        )
        output_name = caller_config[VariantCallingKeys.INDEL_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
            "-V",
            input_name,
            "--select-type-to-include",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def _create_get_other_variants_command(
        cls, caller_config: Dict, library_paths: LibraryPaths
    ) -> List:
        input_name = convert_to_absolute_path(
            join_paths(
                caller_config[VariantCallingKeys.OUTPUT_DIR],
                caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT],
            )
        )
        output_name = caller_config[VariantCallingKeys.OTHER_VARIANTS_OUTPUT]

        command = [
            "gatk",
            "SelectVariants",
            "-R",
            library_paths.REF_FASTA,
            "-V",
            input_name,
            "--select-type-to-exclude",
            "SNP",
            "--select-type-to-exclude",
            "INDEL",
            "-O",
            output_name,
        ]

        return command

    @classmethod
    def call_variants(cls, caller_config: Dict, device: str = "cpu"):
        library_paths = LibraryPaths()
        workdir = caller_config[VariantCallingKeys.OUTPUT_DIR]

        bed_file = (
            caller_config[VariantCallingKeys.BED_FILE]
            if VariantCallingKeys.BED_FILE in caller_config.keys()
            else None
        )

        if device == "cpu":

            splitted_configs = ScatterGather.split_variantcaller_configs(
                caller_config, bed_file=bed_file
            )

            with MemoryHandler() as memory_handler:
                scattered_commands = [
                    cls._create_run_command(
                        caller_config=cfg,
                        library_paths=library_paths,
                        memory_handler=memory_handler,
                    )
                    for cfg in splitted_configs
                ]

                results = run_command_parallel(scattered_commands, cwd=workdir)

                # Check if any of the commands failed
                # If respective region of bam file is empty, results will be None, ignore it
                if any(
                    (result.returncode != 0 and result is not None)
                    for result in results
                ):
                    raise Exception("HaplotypeCaller failed")

            ScatterGather.gather_vcfs(
                splitted_configs,
                output_path=(
                    caller_config[VariantCallingKeys.GVCF_OUTPUT]
                    if caller_config[VariantCallingKeys.OUTPUT_TYPE].lower() == "gvcf"
                    else caller_config[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT]
                ),
                mode=caller_config[VariantCallingKeys.OUTPUT_TYPE],
                cwd=workdir,
            )

            ScatterGather.clean_temp_files(caller_config[VariantCallingKeys.OUTPUT_DIR])

        elif device == "gpu":
            command = cls._create_parabricks_haplotypecaller_command(
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

        if caller_config[VariantCallingKeys.OUTPUT_TYPE] == "VCF":

            output_dir = Path(workdir).absolute()

            cnnscorevariants_command = cls._create_cnnscorevariants_command(
                caller_config=caller_config, library_paths=library_paths
            )
            create_vcfs_index_command = cls._create_vcf_index_command(
                caller_config=caller_config, library_paths=library_paths
            )
            filter_variants_command = cls._create_filter_variants_command(
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

            docker_runner = DockerRunner(device=device)
            docker_runner.run(
                DockerImages.GATK,
                " ".join(cnnscorevariants_command),
                workdir=str(output_dir.parent.parent),
            )                

            run(create_vcfs_index_command, cwd=workdir)
            run(filter_variants_command, cwd=workdir)
            run(get_snp_command, cwd=workdir)
            run(get_indel_command, cwd=workdir)
            run(get_other_variants_command, cwd=workdir)
