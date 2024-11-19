from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict


class _IPipelineStep(ABC):
    log_file: str = None

    @abstractmethod
    def get_output(self) -> str:
        pass

    @abstractmethod
    def get_config(self) -> Dict:
        pass


class _PipelineStep:
    def _get_name_from_path(self, path) -> str:
        return Path(path).name.split(".")[0]
