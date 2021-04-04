from ._formats import FileFormats
from .pipeline_config import PipelineKeys, VariantCallingKeys
from ._pipeline_steps import _PipelineStep
from dataclasses import dataclass
from copy import copy


@dataclass
class VariantCaller(_PipelineStep):
    library: str
    params: Dict
    germline: str = None
    tumor: str = None
    germline_filename: str = None
    tumor_filename: str = None

    def __post_init__(self):
        if self.germline is None and self.germline_filename is None:
            raise Exception("Both germline and germline filename cannot be None.")
        if self.tumor is None and self.tumor_filename is None:
            raise Exception("Both tumor and tumor filename cannot be None.")

    def get_output_file_name(self):
        return self.config[PipelineKeys.VARIANT_CALLING][VariantCallingKeys.OUTPUT]

    def set_inputs(self, tumor: str, germline: str):
        # TODO: maybe an error here if either is set to None
        self.germline_filename = germline
        self.tumor_filename = tumor

    @property
    def config(self):
        output_filename = FileFormats.GATK_SNP_OUTPUT.format(identification=self.name)

        vc_config = {
            VariantCallingKeys.LIBRARY: self.library,
            VariantCallingKeys.NORMAL_SRC: self.germline_filename,
            VariantCallingKeys.TUMOR_SRC: self.tumor_filename,
            VariantCallingKeys.PARAMS: self.params,
            VariantCallingKeys.OUTPUT: output_filename,
        }

        config = {PipelineKeys.VARIANT_CALLING: vc_config}
        return config
    