import os
from subprocess import run

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._pipeline_config import MSICallingKeys
from .._utils import join_paths
from ._msi_callers import _MSICaller


class MSISensorPro(_MSICaller):
    @classmethod
    def _create_msisensorpro_command(
        cls, caller_config: dict, library_paths: LibraryPaths, app_config: AppConfig
    ) -> list:

        normal_input = caller_config[MSICallingKeys.NORMAL_INPUT]
        tumor_input = caller_config[MSICallingKeys.TUMOR_INPUT]
        output = caller_config[MSICallingKeys.OUTPUT]

        command = [
            "msisensor-pro",
            "msi",
            "-d",
            library_paths.MSISENSOR_MICROSATELLITES,
            "-n",
            normal_input,
            "-t",
            tumor_input,
            "-o",
            output,
        ]
        return command

    @classmethod
    def call(cls, caller_config: dict):
        library_paths = LibraryPaths()
        app_config = AppConfig()
        command = cls._create_msisensorpro_command(
            caller_config, library_paths, app_config
        )
        run(command)
