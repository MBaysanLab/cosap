from abc import ABC, abstractmethod

from .._config import AppConfig


class _MSICaller(ABC):
    @abstractmethod
    def call(self):
        pass
