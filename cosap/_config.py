from __future__ import annotations

import multiprocessing
import os
from dataclasses import dataclass
from distutils.util import strtobool
from threading import Lock

import psutil

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
    SNAKEFILE_PATH: str = join_paths(COSAP_PATH, "snakemake_workflows", "Snakefile")
    RAMDISK_PATH: str = os.getenv("COSAP_RAMDISK_PATH", "/dev/shm")

    MAX_THREADS_PER_JOB: int = int(os.getenv("COSAP_THREADS_PER_JOB", multiprocessing.cpu_count()))
    if MAX_THREADS_PER_JOB > multiprocessing.cpu_count():
        MAX_THREADS_PER_JOB = multiprocessing.cpu_count()

    # Max memory is (avaliable memory / number of parallel jobs) in bytes
    MAX_MEMORY_PER_JOBS: int = psutil.virtual_memory().total / (multiprocessing.cpu_count() // MAX_THREADS_PER_JOB)
    WORKDIR: str = os.getcwd()
    
    #Set this to True if you are running cosap on a slurm cluster.
    SLURM_CLUSTER = strtobool(os.getenv("COSAP_SLURM_CLUSTER_MODE", "False"))
    IN_MEMORY_MODE = strtobool(os.getenv("COSAP_IN_MEMORY_MODE", "False"))
    DEVICE = os.getenv("COSAP_DEVICE", "cpu")