import os
from dataclasses import dataclass
from typing import Dict

from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class FastqReader(_IPipelineStep, _PipelineStep):
    filename: str
    read: int
    platform: str = "illumina"

    def __post_init__(self):
        self.read = str(self.read)
        self.filename = os.path.normpath(self.filename)

    def get_output(self) -> str:
        return self.filename

    def get_config(self) -> Dict:
        # TODO: return file details here
        return {}


@dataclass
class BamReader(_IPipelineStep, _PipelineStep):
    filename: str

    def __post_init__(self):
        self.filename = os.path.normpath(self.filename)

    def get_output(self) -> str:
        return self.filename

    def get_config(self) -> Dict:
        # TODO: return file details here
        return {}
