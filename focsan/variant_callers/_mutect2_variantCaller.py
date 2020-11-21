import os
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from ._variantCallers import _Callable, _VariantCaller



class Mutect2VariantCaller(_Callable, _VariantCaller):

    
    @classmethod
    def _create_command(cls, pipeline_config: PipelineConfig):
        pass