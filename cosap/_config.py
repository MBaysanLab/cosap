from __future__ import annotations
import imp

import os
from dataclasses import dataclass
from threading import Lock
from ._utils import join_paths

class _AppConfigMeta(type):
    _instances = {}
    _lock: Lock = Lock()

    def __call__(cls) -> AppConfig:
        with cls._lock:
            if cls not in cls._instances:
                cls._instances[cls] = super().__call__()
        return cls._instances[cls]


@dataclass
class AppConfig(metaclass=_AppConfigMeta):
    COSAP_PATH: str = os.environ.get("COSAP")
    LIBRARY_PATH: str = os.environ.get("COSAP_LIBRARY_PATH")
    SNAKEFILE_PATH: str = join_paths(COSAP_PATH, "snakemake_workflows","Snakefile")

    THREADS: int = 8 #This is number of threads each job can use and not the all available threads
    WORKDIR: str = os.getcwd()
