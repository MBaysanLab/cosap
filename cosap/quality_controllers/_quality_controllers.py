from abc import ABC, abstractmethod


class _QualityController(ABC):
    @abstractmethod
    def run_qualitycontroller(self):
        pass


class _QualityControllable:
    pass
