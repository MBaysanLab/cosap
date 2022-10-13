import os
from subprocess import run
from typing import List

from ._config import AppConfig
from ._formats import FileFormats, OutputFolders
from ._library_paths import LibraryPaths


def join_paths(path: str, *paths) -> str:
    """Joins and normalizes paths according to os standards"""
    return os.path.normpath(os.path.join(path, *paths))
