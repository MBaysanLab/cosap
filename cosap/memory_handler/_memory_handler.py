import os
import shutil
import tempfile
from pathlib import Path

from .._config import AppConfig
from .._utils import join_paths


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

    def get_output_path(self, path: str, save_on_close: bool = True) -> str:
        if not self.in_memory_active:
            return path

        tmp_path = join_paths(self.temp_dir, os.path.basename(path))
        Path(tmp_path).touch()

        if save_on_close:
            self.saving_paths[tmp_path] = path

        return tmp_path

    def get_bam_path(self, path: str) -> str:
        bam_path = self.get_path(path)
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file {bam_path} not found")
        
        # get index path along with bam path. the index can be either in form of .bai or .bam.bai

        if os.path.exists(bam_path + ".bai"):
            bai_path = bam_path + ".bai"
        elif os.path.exists(os.path.splitext(bam_path)[0] + ".bai"):
            bai_path = os.path.splitext(bam_path)[0] + ".bai"
        else:
            raise FileNotFoundError(f"Index file for {bam_path} not found")

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

    def save_tempfiles_to_into_disk(self):
        for tmp_path in self.saving_paths:
            shutil.copy(tmp_path, self.saving_paths[tmp_path])

        for dir in self.opened_dirs:
            dir.cleanup()

    def close_without_saving(self):
        for dir in self.opened_dirs:
            dir.cleanup()
