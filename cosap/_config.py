from __future__ import annotations

import multiprocessing
import os
from dataclasses import dataclass
from distutils.util import strtobool
from threading import Lock

import importlib_resources
from contextlib import ExitStack
import atexit

import psutil


def get_snakefile_path() -> str:
    try:
        file_manager = ExitStack()
        atexit.register(file_manager.close)

        ref = importlib_resources.files("cosap") / "snakemake_workflows/Snakefile"
        path = file_manager.enter_context(importlib_resources.as_file(ref))

        return path
        
    except ModuleNotFoundError:
        # TODO: Handle this more gracefully
        return None


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
    HOME_PATH: str = os.getenv("HOME")

    LIBRARY_PATH: str = os.getenv(
        "COSAP_LIBRARY_PATH",
        os.path.join(HOME_PATH, "cosap_data")
    )
    if not os.path.exists(LIBRARY_PATH):
        os.mkdir(LIBRARY_PATH)

    SNAKEFILE_PATH: str = get_snakefile_path()
    
    RAMDISK_PATH: str = os.getenv("COSAP_RAMDISK_PATH", "/dev/shm")

    MAX_THREADS_PER_JOB: int = int(
        os.getenv("COSAP_THREADS_PER_JOB", multiprocessing.cpu_count())
    )
    if (MAX_THREADS_PER_JOB == 0 or
        MAX_THREADS_PER_JOB > multiprocessing.cpu_count()):
        MAX_THREADS_PER_JOB = multiprocessing.cpu_count()

    # Max memory is (avaliable memory / number of parallel jobs) in bytes
    MAX_MEMORY_PER_JOBS: int = psutil.virtual_memory().total / (
        multiprocessing.cpu_count() // MAX_THREADS_PER_JOB
    )

    WORKDIR: str = os.getcwd()

    # Set this to True if you are running COSAP on a Slurm cluster.
    SLURM_CLUSTER = strtobool(os.getenv("COSAP_SLURM_CLUSTER_MODE", "False"))

    IN_MEMORY_MODE = strtobool(os.getenv("COSAP_IN_MEMORY_MODE", "False"))

    DEVICE = os.getenv("COSAP_DEVICE", "cpu")
