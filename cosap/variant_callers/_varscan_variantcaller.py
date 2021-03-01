import glob
import os
from pathlib import Path
from subprocess import run
from typing import Dict, List, Union

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._variantcallers import _Callable, _VariantCaller


class VarScanVariantCaller(_Callable, _VariantCaller):
    @classmethod
    def _create_samtools_mpileup_command(
        cls, caller_config=Dict, library_paths=LibraryPaths
    ) -> List:
        bam_paths = cls._get_bam_paths(caller_config)

        tumor_sample_name = cls._get_sample_name(bam_paths["tumor_bam_path"])

        snp_output_name = cls._create_output_filename(
            caller_config, sample_name=os.path.join("SNP", tumor_sample_name)
        )
        indel_output_name = cls._create_output_filename(
            caller_config, sample_name=os.path.join("INDEL", tumor_sample_name)
        )

        command = [
            "samtools",
            "mpileup",
            library_paths.REF_DIR,
            "-q",
            "1",
            "-B",
            bam_paths["germline_bam_path"],
            " ",
            bam_paths["tumor_bam_path"],
            "java",
            "-jar",
            library_paths.VARSCAN,
            "somatic",
            "--output-snp",
            snp_output_name,
            "--output-indel",
            indel_output_name,
            "--mpileup",
            "1",
            "--strand-filter",
            "0",
        ]
        return command

    @classmethod
    def _create_process_somatic_command(
        cls,
        caller_config=Dict,
        library_paths=LibraryPaths,
        mpileup_object=Union[str, Path],
    ) -> List:
        command = [
            "java",
            "-jar",
            library_paths.VARSCAN,
            "processSomatic",
            mpileup_object,
            "--min-tumor-freq",
            "0.1",
            "--max-tumor-freq",
            "0.05",
            "--p-value",
            "0.07",
        ]
        return command

    @classmethod
    def call_variants(cls, caller_config: Dict):
        library_paths = LibraryPaths()

        samtools_mpileup = cls._create_samtools_mpileup_command(
            caller_config=caller_config, library_paths=library_paths
        )

        run(samtools_mpileup, cwd=caller_config.VCF_OUTPUT_DIR)

        samtools_pileups = glob.glob(f"*{caller_config.VCF_OUTPUT_DIR}*vcf*")

        for vcf_file in samtools_pileups:
            process_somatic_command = cls._create_process_somatic_command(
                caller_config=caller_config,
                library_paths=library_paths,
                mpileup_object=vcf_file,
            )
            run(process_somatic_command, cwd=caller_config.VCF_OUTPUT_DIR)