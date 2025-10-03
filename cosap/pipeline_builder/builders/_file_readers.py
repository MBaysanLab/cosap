import os
from dataclasses import dataclass
from typing import Dict

import pysam

from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class FastqReader(_IPipelineStep, _PipelineStep):
    filename: str
    read: int
    platform: str = "illumina"
    name: str = None
    next_step: _PipelineStep = None

    def __post_init__(self):
        self.read = str(self.read)
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        if self.name == None:
            self.name = self._get_name_from_path(self.filename)

        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"The file {self.filename} does not exist.")

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
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        if self.name == None:
            self.name = self._get_name_from_path(self.filename)

        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"The file {self.filename} does not exist.")

    def get_output(self) -> str:
        return self.filename

    def get_config(self) -> Dict:
        # TODO: return file details here
        return {}


@dataclass
class VCFReader(_IPipelineStep, _PipelineStep):
    filename: str
    name: str = None

    def __post_init__(self):
        self.filename = os.path.abspath(os.path.normpath(self.filename))

        if self.name == None:
            self.name = self._get_name_from_path(self.filename)

        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"The file {self.filename} does not exist.")

    def get_sample_name(self) -> str:
        return pysam.VariantFile(self.filename).header.samples[0]

    def get_output(self) -> str:
        return self.filename

    def get_config(self) -> Dict:
        # TODO: return file details here
        return {}
