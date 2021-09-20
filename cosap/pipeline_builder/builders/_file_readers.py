import os
from dataclasses import dataclass
from posixpath import lexists
from typing import Dict

from ..._formats import FileFormats, Platforms
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class FastqReader(_IPipelineStep, _PipelineStep):
    filename: str
    read: int
    name: str = None
    platform: str = Platforms.ILLUMINA

    def __post_init__(self):
        self.read = str(self.read)
        self.filename = os.path.normpath(self.filename)
        if self.name is None:
            self.name = f"{os.path.basename(self.filename)}_{self.read}"

    def get_output(self) -> str:
        return self.filename

    def get_config(self) -> Dict:
        # TODO: return file details here
        return {}


@dataclass
class BamReader(_IPipelineStep, _PipelineStep):
    filename: str
    name: str = None

    def __post_init__(self):
        self.filename = os.path.normpath(self.filename)
        if self.name is None:
            self.name = os.path.basename(self.filename)

    def get_output(self) -> str:
        return self.filename

    def get_config(self) -> Dict:
        # TODO: return file details here
        return {}
