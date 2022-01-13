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
        "/media/mae/Elements/arif/hg38_bundle/"
    )
    SNAKEFILE_PATH: str = os.path.normpath(
        "/home/mae/Desktop/cosap/snakemake_workflows/Snakefile"
    )

    THREADS: int = 8 #This is number of threads each job can use and not the all available threads
    WORKDIR: str = os.getcwd()
