from abc import ABC, abstractmethod


class _MitochondrialCaller(ABC):
    @classmethod
    @abstractmethod
    def call(cls, caller_config: dict) -> None:
        pass
