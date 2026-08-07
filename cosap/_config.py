from __future__ import annotations

import atexit
import importlib.resources as importlib_resources
import multiprocessing
import os
import re
from contextlib import ExitStack
from dataclasses import dataclass
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


def _strtobool(value: str) -> bool:
    value = value.strip().lower()
    if value in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif value in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"Invalid truth value: {value}")


def _safe_positive_int_from_env(env_var: str, default: int) -> int:
    raw_value = os.getenv(env_var)
    if raw_value is None:
        return default

    try:
        parsed = int(raw_value)
    except ValueError:
        return default

    return parsed if parsed > 0 else default


def _safe_samtools_memory_from_env(default: str = "768M") -> str:
    value = os.getenv("COSAP_SAMTOOLS_SORT_MEM_PER_THREAD", default)
    if re.match(r"^\d+[KMG]$", value):
        return value
    return default


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
    # Load the .env file
    load_dotenv(dotenv_path=get_cosap_dotenv())

    HOME_PATH = os.getenv("HOME")
    LIBRARY_PATH: str = os.getenv(
        "COSAP_LIBRARY_PATH", os.path.join(HOME_PATH, "cosap_data")
    )
    if not os.path.exists(LIBRARY_PATH):
        os.mkdir(LIBRARY_PATH)

    SNAKEFILE_PATH: str = get_snakefile_path()

    RAMDISK_PATH: str = os.getenv("COSAP_RAMDISK_PATH", "/dev/shm")
    REF_VERSION = os.getenv("COSAP_REF_VERSION", "hg38")

    MAX_THREADS_PER_JOB: int = int(
        os.getenv("COSAP_THREADS_PER_JOB", multiprocessing.cpu_count())
    )
    if MAX_THREADS_PER_JOB == 0 or MAX_THREADS_PER_JOB > multiprocessing.cpu_count():
        MAX_THREADS_PER_JOB = multiprocessing.cpu_count()

    SAMTOOLS_SORT_THREADS: int = _safe_positive_int_from_env(
        "COSAP_SAMTOOLS_SORT_THREADS",
        min(MAX_THREADS_PER_JOB, 3),
    )
    SAMTOOLS_SORT_THREADS = min(SAMTOOLS_SORT_THREADS, MAX_THREADS_PER_JOB)
    SAMTOOLS_SORT_MEM_PER_THREAD: str = _safe_samtools_memory_from_env()

    # Max memory is (avaliable memory / number of parallel jobs) in bytes
    MAX_MEMORY_PER_JOB: int = (
        psutil.virtual_memory().total
        / (multiprocessing.cpu_count() // MAX_THREADS_PER_JOB)
        * 0.8
    )

    WORKDIR: str = os.getcwd()

    # Set this to True if you are running COSAP on a Slurm cluster.
    # SLURM_CLUSTER = strtobool(os.getenv("COSAP_SLURM_CLUSTER_MODE", "False"))

    IN_MEMORY_MODE = _strtobool(os.getenv("COSAP_IN_MEMORY_MODE", "False"))

    DEVICE = os.getenv("COSAP_DEVICE", "cpu")
    VEP_CACHE_VERSION = os.getenv("COSAP_VEP_CACHE_VERSION", "115")
    CONTAINER_RUNTIME = os.getenv("COSAP_CONTAINER_RUNTIME", "docker")
    CONTAINER_ENGINE_BIN = os.getenv("COSAP_APPTAINER_BIN", "apptainer")
    EXECUTOR_URL = os.getenv("COSAP_EXECUTOR_URL", "http://tool-executor:8080").rstrip("/")
    EXECUTOR_TOKEN = os.getenv("COSAP_EXECUTOR_TOKEN", "")
    EXECUTOR_CONNECT_TIMEOUT = float(os.getenv("COSAP_EXECUTOR_CONNECT_TIMEOUT", "5"))
    EXECUTOR_READ_TIMEOUT = float(os.getenv("COSAP_EXECUTOR_READ_TIMEOUT", "30"))
    EXECUTOR_LEASE_SECONDS = int(os.getenv("COSAP_EXECUTOR_LEASE_SECONDS", "90"))
    EXECUTOR_POLL_SECONDS = float(os.getenv("COSAP_EXECUTOR_POLL_SECONDS", "2"))
    EXECUTOR_RETRY_SECONDS = float(os.getenv("COSAP_EXECUTOR_RETRY_SECONDS", "5"))
    EXECUTOR_MAX_WAIT_SECONDS = float(
        os.getenv("COSAP_EXECUTOR_MAX_WAIT_SECONDS", "300")
    )
    EXECUTOR_MAX_RUN_SECONDS = float(
        os.getenv("COSAP_EXECUTOR_MAX_RUN_SECONDS", "86400")
    )
