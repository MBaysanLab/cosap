from abc import ABC, abstractmethod


class _GeneFusionCaller(ABC):
    @abstractmethod
    def call(self):
        pass
