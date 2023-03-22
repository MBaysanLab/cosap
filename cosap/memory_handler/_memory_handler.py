import os
import shutil
import tempfile
from .._utils import join_paths
from .._config import AppConfig
from pathlib import Path


class MemoryHandler:
    def __init__(self):
        self.in_memory_active = AppConfig.IN_MEMORY_MODE
        self.opened_dirs = []
        self.saving_paths = {}
        self.temp_dir = self.get_temp_dir()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            self.close_without_saving()
        else:
            self.save_tempfiles_to_into_disk()

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

    def get_output_path(self, path: str, save_on_close:bool = True) -> str:
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
        _ = self.get_path(path.replace(".bam", ".bai"))

        return bam_path

    def get_temp_dir(self) -> str:
        """
        If in_memory mode is active, create a temporary directory and return the path,
        if not, return the original path.
        """
        if not self.in_memory_active:
            tmp_dir = tempfile.TemporaryDirectory()

        tmp_dir = tempfile.TemporaryDirectory(dir=AppConfig.RAMDISK_PATH)
        self.opened_dirs.append(tmp_dir)

        return tmp_dir.name

    def save_tempfiles_to_into_disk(self):
        for tmp_path in self.saving_paths:
            shutil.copy(tmp_path, self.saving_paths[tmp_path])

        for dir in self.opened_dirs:
            dir.cleanup()

    def close_without_saving(self):
        for dir in self.opened_dirs:
            dir.cleanup()
