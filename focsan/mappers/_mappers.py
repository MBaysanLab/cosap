from abc import ABC, abstractmethod


class _Mapper(ABC):
    @abstractmethod
    def map(self):
        pass


class BWAMapper(_Mapper):
    def map(self):
        pass


class Bowtie2Mapper(_Mapper):
    def map(self):
        pass

class NovoalignMapper(_Mapper):
    def map(self):
        pass

