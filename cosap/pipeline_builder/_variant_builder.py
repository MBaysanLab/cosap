from ._formats import FileFormats
from .pipeline_config import PipelineKeys, VariantCallingKeys
from ._pipeline_steps import _PipelineStep


class VariantCaller(_PipelineStep):
    def __init__(
        self,
        library: str,
        normal_filename: str,
        tumor_filename: str,
        library_params: Dict,
        name: str = None,
    ):
        self.name = name
        config = {
            VariantCallingKeys.LIBRARY: library,
            VariantCallingKeys.NORMAL_SRC: normal_filename,
            VariantCallingKeys.TUMOR_SRC: tumor_filename,
            VariantCallingKeys.PARAMS: library_params,
        }

        self.steps = {PipelineKeys.VARIANT_CALLING: config}
