from abc import ABC, abstractmethod


class _Annotator(ABC):
    @abstractmethod
    def annotate(self):
        pass


class _Annotatable:
    pass
