from abc import ABC, abstractmethod
from typing import Dict
from uuid import uuid4


class _IPipelineStep(ABC):
    @abstractmethod
    def get_output(self) -> str:
        pass

    @abstractmethod
    def get_config(self) -> Dict:
        pass


class _PipelineStep:
    def _get_name(self) -> str:
        return uuid4().hex[:4].upper()
