import os
import shutil
import uuid

from .._config import AppConfig
from .._library_paths import LibraryPaths
from .._utils import join_paths


class MemoryHandler:
    def __init__(self, path_to_save_on_success):
        self.dir_on_ramdisk = join_paths(AppConfig.RAMDISK_PATH, str(uuid.uuid1()))
        self.in_memory_active = AppConfig.IN_MEMORY_MODE
        self.path_to_save_on_success = path_to_save_on_success

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            self.remove_without_saving()
        else:
            self.save_into_drive()

    def get_path(self, path: str, load:bool=True) -> str:
        """
        If in_memory mode is active, load file into ramdisk and return the path,
        if not, return the original path.
        """

        if not self.in_memory_active:
            return path

        dirname = os.path.dirname(path)
        os.makedirs(join_paths(self.dir_on_ramdisk, dirname), exist_ok=True)
        file_path_on_ramdisk = join_paths(self.dir_on_ramdisk, path)

        if os.path.exists(file_path_on_ramdisk):
            return file_path_on_ramdisk

        shutil.copy(path, file_path_on_ramdisk)
        return file_path_on_ramdisk

    def load_tmp_file(self, path):
        if not self.in_memory_active:
            return path

        filename = os.path.basename(path)
        file_path_on_ramdisk = join_paths(AppConfig.RAMDISK_PATH, filename)

        if os.path.exists(file_path_on_ramdisk):
            return file_path_on_ramdisk

        return file_path_on_ramdisk

    def load_ref_into_mem(self):
        library_paths = LibraryPaths()
        ref_path = library_paths.REF_FASTA
        ref_id_path = f"{ref_path}.fai"
        ref_dict_path = ref_path.replace("fasta", "dict")

        ref_path_on_mem = self.load_tmp_file(ref_path)
        _ = self.load_tmp_file(ref_id_path)
        _ = self.load_tmp_file(ref_dict_path)

        return ref_path_on_mem

    def save_into_drive(self):
        if self.in_memory_active:
            shutil.copytree(
                self.dir_on_ramdisk, self.path_to_save_on_success, dirs_exist_ok=True
            )
            shutil.rmtree(self.dir_on_ramdisk)

    def remove_without_saving(self):
        if self.in_memory_active:
            shutil.rmtree(self.dir_on_ramdisk)
