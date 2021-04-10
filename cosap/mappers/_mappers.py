from abc import ABC, abstractmethod



class _Mapper(ABC):
    @abstractmethod
    def map(self):
        pass


class _Mappable:
    pass
