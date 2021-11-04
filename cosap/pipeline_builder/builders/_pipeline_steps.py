from abc import ABC, abstractmethod
from typing import Dict
from pathlib import Path

class _IPipelineStep(ABC):
    @abstractmethod
    def get_output(self) -> str:
        pass

    @abstractmethod
    def get_config(self) -> Dict:
        pass


class _PipelineStep:
    def _get_name_from_path(self, path) -> str:
        return Path(path).name.split(".")[0]
