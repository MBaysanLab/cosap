import glob
import gzip
import os
from abc import ABC, abstractmethod
from collections import defaultdict
from parser import Parser
from typing import Dict, List

from .._pipeline_config import PipelineConfig
from .._utils import join_paths


class _Mapper(ABC):
    @abstractmethod
    def map(self):
        pass


class _Mappable:
    pass
