import os
import shutil
import tempfile
from glob import glob
from pathlib import Path
from typing import Literal

from .._config import AppConfig
from .._utils import join_paths, get_bam_index_path


class MemoryHandler:
    def __init__(self):
        self.in_memory_active: Literal[0, 1] = AppConfig.IN_MEMORY_MODE
        self.opened_dirs: list = []
        self.saving_paths: dict = {}
        self.temp_dir: str = self.get_temp_dir()

    def __enter__(self):
        if not self.temp_dir:
            self.temp_dir = self.get_temp_dir()

        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            self.close_without_saving()
        else:
            self.save_tempfiles_to_disk()

    def get_path(self, path: str) -> str:
        """
        If in_memory mode is active, load file into ramdisk and return the path,
        if not, return the original path.
        """
        if not self.in_memory_active:
            return path

        tmp_path = join_paths(self.temp_dir, os.path.basename(path))

        if os.path.exists(tmp_path):
            return tmp_path
        else:
            shutil.copy(path, self.temp_dir)

        return tmp_path

    def get_output_path(self, path: str, save_on_close: bool = True) -> str:
        if not self.in_memory_active:
            return path

        tmp_path = join_paths(self.temp_dir, os.path.basename(path))
        Path(tmp_path).touch()

        if save_on_close:
            self.saving_paths[tmp_path] = path

        return tmp_path

    def get_bam_path(self, path: str):
        """
        load bam file and its index into ramdisk and return the path
        """

        bam_path = self.get_path(path)
        # get index path along with bam path. the index can be either in form of .bai or .bam.bai
        bai_path = get_bam_index_path(path)
        _ = self.get_path(bai_path)
        return bam_path

    def get_temp_dir(self, dir=None) -> str:
        """
        If in_memory mode is active, create a temporary directory and return the path,
        if not, return the original path.
        """
        if not self.in_memory_active:
            tmp_dir = tempfile.TemporaryDirectory(dir=dir)
            # give access to all users
            os.chmod(tmp_dir.name, 0o777)
        else:
            tmp_dir = tempfile.TemporaryDirectory(dir=AppConfig.RAMDISK_PATH)

        self.opened_dirs.append(tmp_dir)

        return os.path.normpath(tmp_dir.name)

    def save_tempfiles_to_disk(self):
        for tmp_path in self.saving_paths.copy():
            shutil.copy(tmp_path, self.saving_paths.pop(tmp_path))

        for dir in self.opened_dirs.copy():
            dir.cleanup()
            self.opened_dirs.remove(dir)

        self.temp_dir = None

    def close_without_saving(self):
        self.saving_paths = {}

        for dir in self.opened_dirs.copy():
            dir.cleanup()
            self.opened_dirs.remove(dir)
        
        self.temp_dir = None
