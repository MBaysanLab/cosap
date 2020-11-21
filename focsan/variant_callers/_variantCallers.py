import glob
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
    def _get_sample_name(cls, sample_path: str) -> str:

        # This command is not taken from original pipeline code and based on https://github.com/IARCbioinfo/BAM-tricks#extract-sample-name
        command = [
            "samtools",
            "view",
            "-H",
            sample_path,
            "|",
            "grep",
            "'^@RG'",
            "|",
            "sed",
            "'s/.*SM:\([^\t]*\).*/\1/g'",
            "|",
            "uniq",
        ]

        sample_name = run(command, capture_output=True).stdout

    @classmethod
    def _create_output_filename(
        cls, file_info: Dict, pipeline_config: PipelineConfig
    ) -> str:
        return f"{pipeline_config.MAPPER_TYPE}_{pipeline_config.VARIANT_CALLER_TYPE}_{file_info["sample_name"]}"
