from __future__ import annotations

import os
from dataclasses import dataclass
from threading import Lock


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
    LIBRARY_PATH: str = os.path.normpath(
        "/media/monkey/bigmac/omer_research/data/ref/"
    )
    SNAKEFILE_PATH: str = os.path.normpath(
        "/home/monkey/workspace/cosap/snakemake_workflows/Snakefile"
    )
    THREADS: int = 8
    WORKDIR: str = os.getcwd()
