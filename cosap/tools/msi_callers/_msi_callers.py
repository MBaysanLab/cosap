from abc import ABC, abstractmethod


class _MSICaller(ABC):
    @abstractmethod
    def call(self):
        pass
