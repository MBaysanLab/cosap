import os
import requests
import hashlib
from enum import Enum
import functools

from .._utils import join_paths
from .._config import AppConfig


def _get_file_hash_md5(file_path: str, block_size: int = 2**20) -> str:
    hasher = hashlib.md5()
    file = open(file_path, "rb")

    block = file.read(block_size)
    while block:
        hasher.update(block)
        block = file.read(block_size)
    
    return hasher.hexdigest()


def _reporthook(downloaded_size, total_file_size):
    downloaded_ratio = downloaded_size / total_file_size
    downloaded_percentage = downloaded_ratio * 100 if downloaded_ratio < 1 else 100.0
    print(f"\r  {downloaded_percentage:.2f}%", end="")
    if downloaded_percentage == 100.0:
        print(end="\n")


class DownloadableFileStatus(Enum):
    TO_BE_DOWNLOADED = "to_be_downloaded"
    IGNORED = "ignored"

    PRESENT_FULLY = "present_fully"
    PRESENT_PARTIALLY = "present_partially"

    DOWNLOAD_FAILED = "download_failed"


class DownloadableFile:
    def __init__(
        self,
        url: str,
        filename: str = None,
        size: int = None,
        md5: str = None,
        requiring_steps: set = None
    ):
        self.url = url

        self.filename = filename if filename else url.split("/")[-1]

        self.expected_size = int(size) if size else int(requests.get(url).headers.get("content-length"))

        self.md5 = md5

        self.download_path = join_paths(AppConfig.LIBRARY_PATH, self.filename)
        
        self.status = DownloadableFileStatus.TO_BE_DOWNLOADED

        self.downloaded_size = 0
        if os.path.isfile(self.download_path):
            self.downloaded_size = os.path.getsize(self.download_path)

        self.requiring_steps = set(requiring_steps) if requiring_steps else set()

    
    def ignored(self):
        return self.status == DownloadableFileStatus.IGNORED
    

    def ignorable(inner_function):
        @functools.wraps(inner_function)
        def check_ignore(self, *args, **kwargs):

            if self.ignored():
                return
            
            result = inner_function(self, *args, **kwargs)
            return result
        
        return check_ignore


    @ignorable
    def check_file_integrity(self):
        print(f"{self.filename}: ", end="", flush=True)
        if os.path.isfile(self.download_path):
            self.downloaded_size = os.path.getsize(self.download_path)

            print("File exists")
            if self.check_size(delete_on_too_large=True):
                if self.check_hash(delete_on_mismatch=True):
                    self.status = DownloadableFileStatus.PRESENT_FULLY
                else:
                    self.status = DownloadableFileStatus.TO_BE_DOWNLOADED
        else:
            print("To be downloaded")

    
    @ignorable
    def download(self, retry_count=0):
        if self.status == DownloadableFileStatus.PRESENT_FULLY:
            return
        
        print(f"{self.filename}: Downloading ({self.expected_size:,} bytes)")

        chunk_size = 8192
        
        headers = {"Range": f"bytes={self.downloaded_size}-"}

        if self.downloaded_size < self.expected_size:
            with requests.get(self.url, stream=True, headers=headers) as response:
                with open(self.download_path, "ab") as output_file:
                    for i, chunk in enumerate(response.iter_content(chunk_size=chunk_size)):
                        output_file.write(chunk)
                        _reporthook(self.downloaded_size + ((i+1)*chunk_size), self.expected_size)

        if self.check_size(delete_on_too_large=True) and self.check_hash(delete_on_mismatch=True):
            self.status = DownloadableFileStatus.PRESENT_FULLY
        else:
            if retry_count >= 2:
                print("  This file cannot be downloaded at this time.")
                print("  Please try again in a few hours.")
                self.status = DownloadableFileStatus.DOWNLOAD_FAILED
            else:
                print("Retrying...")
                self.download(retry_count=retry_count+1)
    

    @ignorable
    def check_size(self, delete_on_too_large=False):
        print("  Checking size: ", end="", flush=True)
        self.downloaded_size = os.path.getsize(self.download_path)

        if not self.expected_size:
            print("UNSPECIFIED")
            return False

        if self.downloaded_size < self.expected_size:
            self.status = DownloadableFileStatus.PRESENT_PARTIALLY
            print("FAILED")
            return False

        if self.downloaded_size > self.expected_size:
            print("FAILED")
            if delete_on_too_large:
                print(f"File larger than expected. Deleting {self.filename}.")
                self.delete_file()
            return False

        print("OK")
        return True
        

    @ignorable
    def check_hash(self, delete_on_mismatch=False):
        print("  Checking hash: ", end="", flush=True)

        if not self.md5:
            print("UNSPECIFIED")
            return True
        
        if self.md5 != _get_file_hash_md5(self.download_path):
            print("FAILED")
            if delete_on_mismatch:
                print(f"Deleting {self.filename}.")
                self.delete_file()
            return False

        print("OK")
        return True
    

    @ignorable
    def delete_file(self):
        if os.path.isfile(self.download_path):
            os.remove(self.download_path)
        
        self.downloaded_size = 0
        self.status = DownloadableFileStatus.TO_BE_DOWNLOADED


    def is_required_for(self, steps) -> bool:
        """
        Returns whether the file is required for the given steps.
        Return true if no steps are specified.
        """
        if len(steps) == 0:
            return
        
        requiring_steps = self.requiring_steps & steps

        return 0 < len(requiring_steps)
    
