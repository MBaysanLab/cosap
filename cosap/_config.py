from __future__ import annotations

import atexit
import importlib.resources as importlib_resources
import multiprocessing
import os
from contextlib import ExitStack
from dataclasses import dataclass
from distutils.util import strtobool
from threading import Lock

import psutil
from dotenv import load_dotenv


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
    
def get_cosap_dotenv():
    """
    Create and return a .env file in ~/.cosap directory if it does not exist.
    """

    dotenv_path = os.path.join(os.path.expanduser("~"), ".cosap", ".env")
    os.makedirs(os.path.dirname(dotenv_path), exist_ok=True)

    if not os.path.exists(dotenv_path):
        with open(dotenv_path, "w") as f:
            f.write("")

    return dotenv_path

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

    # Load the .env file
    load_dotenv(
        dotenv_path=get_cosap_dotenv()
    )

    HOME_PATH = os.getenv("HOME")
    LIBRARY_PATH: str = os.getenv(
        "COSAP_LIBRARY_PATH",
        os.path.join(HOME_PATH, "cosap_data")
    )
    if not os.path.exists(LIBRARY_PATH):
        os.mkdir(LIBRARY_PATH)

    SNAKEFILE_PATH: str = get_snakefile_path()
    
    RAMDISK_PATH: str = os.getenv("COSAP_RAMDISK_PATH", "/dev/shm")
    REF_VERSION = os.getenv("COSAP_REF_VERSION", "hg38")

    MAX_THREADS_PER_JOB: int = int(
        os.getenv("COSAP_THREADS_PER_JOB", multiprocessing.cpu_count())
    )
    if (MAX_THREADS_PER_JOB == 0 or
        MAX_THREADS_PER_JOB > multiprocessing.cpu_count()):
        MAX_THREADS_PER_JOB = multiprocessing.cpu_count()

    # Max memory is (avaliable memory / number of parallel jobs) in bytes
    MAX_MEMORY_PER_JOB: int = psutil.virtual_memory().total / (
        multiprocessing.cpu_count() // MAX_THREADS_PER_JOB
    ) * 0.8

    WORKDIR: str = os.getcwd()

    # Set this to True if you are running COSAP on a Slurm cluster.
    # SLURM_CLUSTER = strtobool(os.getenv("COSAP_SLURM_CLUSTER_MODE", "False"))

    IN_MEMORY_MODE = strtobool(os.getenv("COSAP_IN_MEMORY_MODE", "False"))

    DEVICE = os.getenv("COSAP_DEVICE", "cpu")
