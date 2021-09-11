from abc import ABC, abstractmethod


class _Preprocessor(ABC):
    @abstractmethod
    def run_preprocessor(self):
        pass

class _PreProcessable:
    pass
