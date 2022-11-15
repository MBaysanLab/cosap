from abc import ABC, abstractmethod

from .._config import AppConfig


class _Mapper(ABC):
    @classmethod
    def _samtools_sort_command(cls, app_config: AppConfig, output_path: str):
        command = [
            "samtools",
            "sort",
            "-@",
            str(app_config.THREADS),
            "-o",
            output_path,
            "-",
        ]
        return command

    @classmethod
    def _samtools_index_command(cls, app_config: AppConfig, input_path: str):
        command = [
            "samtools",
            "index",
            input_path,
            "-@",
            str(app_config.THREADS),
        ]
        return command

    @abstractmethod
    def map(self):
        pass


class _Mappable:
    pass
