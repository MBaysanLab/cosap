from copy import copy
from dataclasses import dataclass
from typing import Dict

from ..._formats import FileFormats
from ..._pipeline_config import PipelineKeys, VariantCallingKeys, MappingKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class VariantCaller(_IPipelineStep, _PipelineStep):
    library: str
    germline: str
    tumor: str
    params: Dict
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = self._get_name()
        if VariantCallingKeys.GERMLINE_SAMPLE_NAME not in self.params:
            self.params[VariantCallingKeys.GERMLINE_SAMPLE_NAME] = "normal_sample"
        if VariantCallingKeys.TUMOR_SAMPLE_NAME not in self.params:
            self.params[VariantCallingKeys.TUMOR_SAMPLE_NAME] = "tumor_sample"

    def get_output(self):
        config = self.get_config()
        return config[PipelineKeys.VARIANT_CALLING][self.name][
            VariantCallingKeys.SNP_OUTPUT
        ]

    def get_config(self) -> Dict:
        unfiltered_variants_output_filename = FileFormats.GATK_UNFILTERED_OUTPUT.format(
            germline_identification=self.germline.name,
            tumor_identification=self.tumor.name,
            algorithm=self.library,
        )
        snp_output_filename = FileFormats.GATK_SNP_OUTPUT.format(
            germline_identification=self.germline.name,
            tumor_identification=self.tumor.name,
            algorithm=self.library,
        )
        indel_output_filename = FileFormats.GATK_INDEL_OUTPUT.format(
            germline_identification=self.germline.name,
            tumor_identification=self.tumor.name,
            algorithm=self.library,
        )
        other_variants_output_filename = FileFormats.GATK_OTHER_VARIANTS_OUTPUT.format(
            germline_identification=self.germline.name,
            tumor_identification=self.tumor.name,
            algorithm=self.library,
        )

        vc_config = {
            self.name: {
                VariantCallingKeys.LIBRARY: self.library,
                VariantCallingKeys.GERMLINE_INPUT: self.germline.get_output(),
                VariantCallingKeys.TUMOR_INPUT: self.tumor.get_output(),
                VariantCallingKeys.PARAMS: self.params,
                VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT: unfiltered_variants_output_filename,
                VariantCallingKeys.SNP_OUTPUT: snp_output_filename,
                VariantCallingKeys.INDEL_OUTPUT: indel_output_filename,
                VariantCallingKeys.OTHER_VARIANTS_OUTPUT: other_variants_output_filename,
            }
        }

        config = {PipelineKeys.VARIANT_CALLING: vc_config}
        return config
