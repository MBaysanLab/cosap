from copy import copy
from dataclasses import dataclass
from typing import Dict

from .._formats import FileFormats
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from .._pipeline_config import PipelineKeys, VariantCallingKeys


@dataclass
class VariantCaller(_IPipelineStep, _PipelineStep):
    library: str
    germline: str
    tumor: str
    params: Dict
    name: str = None

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.tumor.name}-{self.germline.name}"

    def get_output(self):
        config = self.get_config()
        return config[PipelineKeys.VARIANT_CALLING][VariantCallingKeys.OUTPUT]

    def get_config(self) -> Dict:
        output_filename = FileFormats.GATK_SNP_OUTPUT.format(identification=self.name)

        vc_config = {
            VariantCallingKeys.LIBRARY: self.library,
            VariantCallingKeys.GERMLINE: self.germline.get_output(),
            VariantCallingKeys.TUMOR: self.tumor.get_output(),
            VariantCallingKeys.PARAMS: self.params,
            VariantCallingKeys.OUTPUT: output_filename,
        }

        config = {PipelineKeys.VARIANT_CALLING: vc_config}
        return config
