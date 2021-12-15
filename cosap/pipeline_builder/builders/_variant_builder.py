from copy import copy
from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats
from ..._pipeline_config import MappingKeys, PipelineKeys, VariantCallingKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class VariantCaller(_IPipelineStep, _PipelineStep):
    library: str
    germline: str
    tumor: str
    params: Dict
    name: str = None
    key: str = PipelineKeys.VARIANT_CALLING

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.germline.name}-{self.tumor.name}_{self.library}"
        if VariantCallingKeys.GERMLINE_SAMPLE_NAME not in self.params:
            self.params[VariantCallingKeys.GERMLINE_SAMPLE_NAME] = "normal_sample"
        if VariantCallingKeys.TUMOR_SAMPLE_NAME not in self.params:
            self.params[VariantCallingKeys.TUMOR_SAMPLE_NAME] = "tumor_sample"

    def get_output(self):
        config = self.get_config()
        return config[self.key][self.name][VariantCallingKeys.SNP_OUTPUT]

    def get_config(self) -> Dict:
        unfiltered_variants_output_filename = FileFormats.GATK_UNFILTERED_OUTPUT.format(
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
            VariantCallingKeys.SNAKEMAKE_OUTPUT: FileFormats.GATK_SNP_OUTPUT.format(
                identification="{identification}"
            ),
            self.name: {
                VariantCallingKeys.LIBRARY: self.library,
                VariantCallingKeys.GERMLINE_INPUT: self.germline.get_output(),
                VariantCallingKeys.TUMOR_INPUT: self.tumor.get_output(),
                VariantCallingKeys.PARAMS: self.params,
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: unfiltered_variants_output_filename,
                VariantCallingKeys.SNP_OUTPUT: snp_output_filename,
                VariantCallingKeys.INDEL_OUTPUT: indel_output_filename,
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: other_variants_output_filename,
            },
        }

        config = {self.key: vc_config}
        return config
