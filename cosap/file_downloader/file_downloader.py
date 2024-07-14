import os
import json
from distutils.util import strtobool

import yaml

from .._config import AppConfig
from .._utils import join_paths, swap_dict_keys_and_values, prompt_continue

from .downloadable_file import DownloadableFile, DownloadableFileStatus
from .files import downloadable_files


class FileDownloader:
    def __init__(self) -> None:
        self.download_directory: str = AppConfig.LIBRARY_PATH
        self.checked_hashes_yaml_path: str = join_paths(self.download_directory, "checked_hashes.yaml")
        
        self.checked_files: dict = self.read_checked_files()

    
    def download_files(self, steps: set = set(), confirm: bool = True):
        """
        Downloads files required for steps specified.
        Downloads all possible files if empty set is passed.
        """

        # TODO: Check available space
        print(f"Files will be downloaded to ` {AppConfig.LIBRARY_PATH} `. The directory will be created if it does not exist.")
        print("To change the download location, set the environment variable ` COSAP_LIBRARY_PATH ` to the desired path.")

        if 0 == len(self.checked_files) and confirm:
            if not prompt_continue():
                print("No files downloaded.")
                return
        
        os.makedirs(AppConfig.LIBRARY_PATH, exist_ok=True)

        files_to_be_downloaded = []
        for i, downloadable_file in enumerate(downloadable_files):

            if downloadable_file.filename in self.checked_files.keys():
                present_hash_value = self.checked_files[downloadable_file.filename]

                if len(present_hash_value) < 16:  # User specified "ignore"
                    downloadable_file.status = DownloadableFileStatus.IGNORED
                elif downloadable_file.md5 == present_hash_value:  # TODO: Delete file if target hash is different
                    downloadable_file.status = DownloadableFileStatus.PRESENT_FULLY

            if (
                len(steps) == 0 or
                downloadable_file.is_required_for(steps)
            ):
                files_to_be_downloaded.append(downloadable_file)

        total_files_count = len(files_to_be_downloaded)

        if total_files_count == 0:
            print("-" * 40)
            print("No new files required.")
            return
        
        incomplete_files = []
        previously_present_files = []
        for file in files_to_be_downloaded:
            if file.status == DownloadableFileStatus.PRESENT_FULLY:
                previously_present_files.append(file)
            elif not file.ignored():
                incomplete_files.append(file)
        
        if incomplete_files:
            print("-" * 40)
            print("Downloading files...")

        for i, file in enumerate(incomplete_files):
            print(f"# File {i+1}/{len(incomplete_files)}: ", end="", flush=True)
            file.download()
            if file.status == DownloadableFileStatus.PRESENT_FULLY:
                self.checked_files[file.filename] = file.md5
        
        statuses = {file.filename: file.status for file in files_to_be_downloaded}
        
        print("-" * 40)
        
        statuses_list = list(statuses.values())
        already_present_count = len(previously_present_files)
        ignored_count = statuses_list.count(DownloadableFileStatus.IGNORED)
        downloaded_count = statuses_list.count(DownloadableFileStatus.PRESENT_FULLY) - already_present_count + ignored_count
        download_blocked_count = statuses_list.count(DownloadableFileStatus.DOWNLOAD_FAILED)

        if already_present_count:
            print(f"{already_present_count}/{total_files_count} files were already present.")

        if downloaded_count:
            print(f"{downloaded_count}/{total_files_count} files were downloaded successfully.")

        if download_blocked_count:
            print(f"{download_blocked_count}/{total_files_count} files could not be downloaded:")
            [print("-", filename) for filename in statuses.keys() if statuses[filename] == DownloadableFileStatus.DOWNLOAD_FAILED]

        if ignored_count:
            print(f"{ignored_count}/{total_files_count} files were ignored.")
        
        missing_or_incomplete_count = total_files_count - (already_present_count + downloaded_count + ignored_count)
        if missing_or_incomplete_count:
            print(f"{missing_or_incomplete_count}/{total_files_count} are missing or incomplete.")
        else:
            print("All files are present.")
        
        self.write_checked_files()


    
    def read_checked_files(self) -> dict:
        """
        Read 'checked_hashes.yaml' present in the COSAP data directory.
        Returns a dictionary in the form {filename: hash}.
        """

        checked_hashes = {}
        if os.path.exists(self.checked_hashes_yaml_path):
            with open(self.checked_hashes_yaml_path) as file:
                checked_hashes = yaml.safe_load(file)

            if checked_hashes is None:
                checked_hashes = {}
        
        # It is more readable to have the hashes line up in the YAML file
        checked_files = swap_dict_keys_and_values(checked_hashes)

        existent_files = os.listdir(self.download_directory)
        for file in checked_files.copy().keys():
            if file not in existent_files:
                checked_files.pop(file)

        # TODO: If no files can be found out of many, ask if directory might be mixed up

        return checked_files
        

    def write_checked_files(self):
        """Write to 'checked_hashes.yaml' present in the COSAP data directory."""
        checked_hashes = swap_dict_keys_and_values(self.checked_files)
        with open(self.checked_hashes_yaml_path, "w") as file:
            yaml.safe_dump(checked_hashes, file)

