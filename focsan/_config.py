from __future__ import annotations
from dataclasses import dataclass
from threading import Lock
import os


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
    LIBRARY_PATH: str = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/".replace("/", os.sep)
