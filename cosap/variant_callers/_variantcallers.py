from abc import ABC, abstractmethod


class _VariantCaller(ABC):
    @abstractmethod
    def call_variants(self):
        pass


class _Callable:
    pass
