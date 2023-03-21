import shutil
import tempfile
import os

from .._config import AppConfig

class MemoryHandler:
    def __init__(self, path_to_save_on_success):
        self.in_memory_active = AppConfig.IN_MEMORY_MODE
        self.path_to_save_on_success = path_to_save_on_success
        self.opened_files = {}
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

    def get_path(
        self, path: str, save_path: str = None
    ) -> str:
        """
        If in_memory mode is active, load file into ramdisk and return the path,
        if not, return the original path.
        """
        if not self.in_memory_active:
            return path

        if path in self.opened_files:
            return self.opened_files[path]
        else:
            tmp_file = tempfile.NamedTemporaryFile(prefix=os.path.basename(path), dir = self.temp_dir)
            self.opened_files[path] = tmp_file.name
            if save_path:
                self.saving_paths[tmp_file.name] = save_path

        shutil.copy(path, tmp_file.name)

        return tmp_file.name

    def get_output_path(self, path: str) -> str:
        if not self.in_memory_active:
            return path

        tmp_file = tempfile.NamedTemporaryFile(prefix=os.path.basename(path), dir = self.temp_dir)
        self.opened_files[path] = tmp_file.name

        return tmp_file.name
    

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

        tmp_dir = tempfile.TemporaryDirectory(dir = AppConfig.RAMDISK_PATH)
        self.opened_dirs.append(tmp_dir)

        return tmp_dir.name
    
    def save_tempfiles_to_into_disk(self):
        for file in self.saving_paths:
            shutil.copy(file, self.saving_paths[file])
        
        for file in self.opened_files.values():
            file.close()
        
        for dir in self.opened_dirs.values():
            dir.cleanup()

    def close_without_saving(self):
        for file in self.opened_files.values():
            file.close()

        for dir in self.opened_dirs.values():
            dir.cleanup()