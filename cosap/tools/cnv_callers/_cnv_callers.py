from abc import ABC, abstractmethod

from ..._config import AppConfig


class _CNVCaller(ABC):
    @abstractmethod
    def call(self):
        pass
