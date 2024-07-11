import os
import requests
import hashlib
import json
from enum import Enum


cosap_library_path = os.getenv(
    "COSAP_LIBRARY_PATH",
    default=os.path.join(os.getenv("HOME"), "cosap_data")
)


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


class DownloadFileStatus(Enum):
    to_be_downloaded = "to_be_downloaded"
    present_fully = "present_fully"
    present_partially = "present_partially"
    downloaded_fully = "downloaded_fully"
    downloaded_partially = "downloaded_partially"
    download_blocked = "download_blocked"


class DownloadFile:
    def __init__(self, url, filename=None, size=None, md5=None):
        self.url = url

        self.filename = filename if filename else url.split("/")[-1]

        self.size = size if size else int(requests.get(url).headers.get("content-length"))

        self.md5 = md5

        self.download_path = os.path.join(cosap_library_path, self.filename)
        
        self.status = DownloadFileStatus.to_be_downloaded
        self.downloaded_size = 0

        print(f"{self.filename}: ", end="", flush=True)
        if os.path.isfile(self.download_path):
            self.downloaded_size = os.path.getsize(self.download_path)

            print("File exists")
            if self.check_size():
                if self.check_hash(delete_on_mismatch=True):
                    self.status = DownloadFileStatus.present_fully
                else:
                    self.status = DownloadFileStatus.to_be_downloaded
            else:
                self.status = DownloadFileStatus.present_partially
        else:
            print("To be downloaded")

    

    def download(self, retry_count=0):
        if self.status in (DownloadFileStatus.present_fully, DownloadFileStatus.downloaded_fully):
            return
        
        print(f"{self.filename}: Downloading ({self.size:,} bytes)")

        chunk_size = 8192
        
        headers = {"Range": f"bytes={self.downloaded_size}-"}

        with requests.get(self.url, stream=True, headers=headers) as response:
            with open(self.download_path, "ab") as output_file:
                for i, chunk in enumerate(response.iter_content(chunk_size=chunk_size)):
                    output_file.write(chunk)
                    _reporthook(self.downloaded_size + ((i+1)*chunk_size), self.size)

        if self.check_size() and self.check_hash(delete_on_mismatch=True):
            self.status = DownloadFileStatus.downloaded_fully
        else:
            if retry_count >= 2:
                print("  This file cannot be downloaded at this time.")
                print("  Please try again in a few hours.")
                self.status = DownloadFileStatus.download_blocked
            else:
                print("Retrying...")
                self.download(retry_count=retry_count+1)
    

    def check_size(self):
        print("  Checking size: ", end="", flush=True)

        if not self.size:
            print("UNSPECIFIED")
            return False

        present_file_size = os.path.getsize(self.download_path)
        if self.size != present_file_size:
            print("FAILED")
            return False

        print("OK")
        return True
        

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
    

    def delete_file(self):
        if os.path.isfile(self.download_path):
            os.remove(self.download_path)
        
        self.downloaded_size = 0
        self.status = DownloadFileStatus.to_be_downloaded


if __name__ == "__main__":
    # TODO: Check available space
    print(f"Files will be downloaded to ` {cosap_library_path} `. The directory will be created if it does not exist.")
    print("To change the download location, set the environment variable ` COSAP_LIBRARY_PATH ` to the desired path.")

    if input("Continue? [Y/n]: ") not in ["", "y", "Y"]:
        print("No files downloaded.")
        exit()
    
    print("-" * 40)
    
    os.makedirs(cosap_library_path, exist_ok=True)

    file_info_list = json.load(open("files.json"))["files"]
    total_files_count = len(file_info_list)

    print("Loading file info...")
    files_to_be_downloaded = []
    for i, file in enumerate(file_info_list):
        print(f"# File {i+1}/{total_files_count}: ", end="", flush=True)
        files_to_be_downloaded.append(DownloadFile(**file))
    
    incomplete_files = [
        file for file in files_to_be_downloaded
        if file.status not in (
            DownloadFileStatus.present_fully,
            DownloadFileStatus.downloaded_fully
        )
    ]

    if incomplete_files:
        print("-" * 40)
        print("Downloading files...")

    for i, file in enumerate(incomplete_files):
        print(f"# File {i+1}/{len(incomplete_files)}: ", end="", flush=True)
        file.download()
    
    statuses = {file.filename: file.status for file in files_to_be_downloaded}
    
    print("-" * 40)
    
    statuses_list = list(statuses.values())
    already_present_count = statuses_list.count(DownloadFileStatus.present_fully)
    downloaded_count = statuses_list.count(DownloadFileStatus.downloaded_fully)
    download_blocked_count = statuses_list.count(DownloadFileStatus.download_blocked)

    if already_present_count:
        print(f"{already_present_count}/{total_files_count} files were already present.")

    if downloaded_count:
        print(f"{downloaded_count}/{total_files_count} files were downloaded successfully.")

    if download_blocked_count:
        print(f"{download_blocked_count}/{total_files_count} files could not be downloaded:")
        [print("-", filename) for filename in statuses.keys() if statuses[filename] == DownloadFileStatus.download_blocked]
    
    if already_present_count + downloaded_count < total_files_count:
        print(f"{total_files_count - (already_present_count + downloaded_count)}/{total_files_count} are missing or incomplete.")
    else:
        print("All files are present.")
