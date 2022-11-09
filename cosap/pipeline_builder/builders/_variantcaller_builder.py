import imp
from copy import copy
from dataclasses import dataclass, field
from subprocess import PIPE, STDOUT, Popen
from typing import Dict

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import (
    DefaultValues,
    MappingKeys,
    PipelineKeys,
    VariantCallingKeys,
)
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class VariantCaller(_IPipelineStep, _PipelineStep):
    PY2_PACKAGES = ["strelka"]

    library: str
    params: dict = field(default_factory=dict)
    germline: str = None
    tumor: str = None
    name: str = None
    key: str = PipelineKeys.VARIANT_CALLING

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.germline.name}-{self.tumor.name}_{self.library}"

        if VariantCallingKeys.GERMLINE_SAMPLE_NAME not in self.params and self.germline:
            self.params[
                VariantCallingKeys.GERMLINE_SAMPLE_NAME
            ] = self._get_sample_name_from_bam(self.germline)

    def _get_sample_name_from_bam(self, bam) -> str:
        cmd = f"samtools view -H {bam} | grep '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/g' | uniq"
        ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        sample_name = str(ps.communicate()[0])
        return str(sample_name)

    def get_output(self):
        config = self.get_config()
        return config[self.key][self.name][VariantCallingKeys.SNP_OUTPUT]

    def get_config(self) -> Dict:
        unfiltered_variants_output_filename = FileFormats.GATK_UNFILTERED_OUTPUT.format(
            identification=self.name
        )
        filtered_variants_output_filename = FileFormats.GATK_FILTERED_OUTPUT.format(
            identification=self.name
        )
        snp_output_filename = FileFormats.GATK_SNP_OUTPUT.format(
            identification=self.name
        )
        indel_output_filename = FileFormats.GATK_INDEL_OUTPUT.format(
            identification=self.name
        )
        other_variants_output_filename = FileFormats.GATK_OTHER_VARIANTS_OUTPUT.format(
            identification=self.name
        )
        vc_config = {
            self.name: {
                VariantCallingKeys.LIBRARY: self.library,
                VariantCallingKeys.PARAMS: self.params,
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING,
                    self.library,
                    unfiltered_variants_output_filename,
                ),
                VariantCallingKeys.FILTERED_VARIANTS_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING,
                    self.library,
                    filtered_variants_output_filename,
                ),
                VariantCallingKeys.SNP_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library, snp_output_filename
                ),
                VariantCallingKeys.INDEL_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library, indel_output_filename
                ),
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: join_paths(
                    OutputFolders.VARIANT_CALLING,
                    self.library,
                    other_variants_output_filename,
                ),
                VariantCallingKeys.OUTPUT_DIR: join_paths(
                    OutputFolders.VARIANT_CALLING, self.library
                ),
            },
        }

        if self.tumor is not None:
            vc_config[self.name][
                VariantCallingKeys.TUMOR_INPUT
            ] = self.tumor.get_output()

        if self.germline is not None:
            vc_config[self.name][
                VariantCallingKeys.GERMLINE_INPUT
            ] = self.germline.get_output()

        config = {self.key: vc_config}
        return config