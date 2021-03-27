import glob
import json
import os
from abc import ABC, abstractmethod
from subprocess import run
from typing import Dict, List

from .._pipeline_config import VariantCallingKeys
from .._utils import join_paths


class _VariantCallable(ABC):
    @abstractmethod
    def call_variants(self):
        pass


class _Callable:
    pass
