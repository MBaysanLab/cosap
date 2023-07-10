from abc import ABC, abstractmethod

from .._config import AppConfig


class _GeneFusionCaller(ABC):
    @abstractmethod
    def call(self):
        pass
