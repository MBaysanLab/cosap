import glob
import json
import os
import tempfile
from abc import ABC, abstractmethod
from subprocess import run
from typing import Dict, List

from .._pipeline_config import PipelineConfig
from .._utils import join_paths


class _VariantCaller(ABC):
    @abstractmethod
    def call_variants(self):
        pass


class _Callable:
    @classmethod
    def _list_vcf_files(cls, file_path: str) -> List:
        return glob.glob(os.path.join(file_path, "*.vcf*"))

    @classmethod
    def _get_bam_paths(cls, pipeline_config: PipelineConfig) -> Dict:

        variant_calling_data = json.loads(pipeline_config.USER_CONFIG)[
            "variant-calling"
        ]["data"]
        file_paths = {
            "germline_bam": variant_calling_data["germline_bam_path"],
            "turmor_bam": variant_calling_data["tumor_bam_path"],
        }
        return file_paths

    @classmethod
    def _get_sample_name(cls, sample_path: str) -> str:

        # This command is not taken from original pipeline code and based on https://github.com/IARCbioinfo/BAM-tricks#extract-sample-name
        command = [
            "samtools",
            "view",
            "-H",
            sample_path,
            "|",
            "grep",
            "^@RG",
            "|",
            "sed",
            r"s/.*SM:\([^\t]*\).*/\1/g",
            "|",
            "uniq",
        ]

        sample_name = run(command, capture_output=True).stdout
        return sample_name

    @classmethod
    def _create_output_filename(
        cls, pipeline_config: PipelineConfig, sample_name: str
    ) -> str:
        return f"{pipeline_config.MAPPER_TYPE}_{pipeline_config.VARIANT_CALLER_TYPE}_{sample_name}.vcf"
