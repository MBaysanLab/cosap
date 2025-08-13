from abc import ABC, abstractmethod


class _CNVCaller(ABC):
    @abstractmethod
    def call(self):
        pass
