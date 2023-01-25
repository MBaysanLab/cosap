import os
from subprocess import run
from typing import List

from .._config import AppConfig
from .._formats import FileFormats, FolderedOutputs
from .._library_paths import LibraryPaths
from .._utils import join_paths
from glob import glob


def merge_bam_files(bam_list: List, output_path: str):
    command = ["samtools", "merge", output_path, *bam_list]
    run(command)


def merge_vcf_files(vcf_list: List, output_path: str):
    command = ["picard", "GatherVcfs", "-I", *vcf_list, "-O", output_path]
    run(command)


def convert_region_file_to_bed(interval_file_path: str):
    output_filename = interval_file_path.split(".")[0]
    command = [
        "picard",
        "IntervalListToBed",
        "-I",
        interval_file_path,
        "-O",
        f"{output_filename}.bed",
    ]
    run(command)

def create_gatk_intervals(
    library_path: LibraryPaths, scatter_count: int, output_dir: str
):  
    command = [
        "gatk",
        "SplitIntervals",
        "-R",
        library_path.REF_FASTA,
        "--scatter-count",
        str(scatter_count),
        "-O",
        output_dir,
    ]
    run(command)

    for interval_file in os.listdir(output_dir):
        convert_region_file_to_bed(join_paths(output_dir,interval_file))

def get_region_file_list(file_type="bed") -> list:
    library_paths = LibraryPaths()
    app_config = AppConfig()

    region_files_dir = FolderedOutputs.REGIONS_FILE_OUTPUT
    if not os.path.exists(region_files_dir):
        create_gatk_intervals(
            library_path=library_paths,
            scatter_count=app_config.MAX_THREADS_PER_JOB,
            output_dir=region_files_dir,
        )
    
    interval_files = glob(f"{region_files_dir}/*.{file_type}")
    return interval_files

def split_bam_by_intervals(bam_path: str) -> list:
    """Splits bam by the number of threads and returns list of names of splitted bams"""

    app_config = AppConfig()
    output_list = []
    threads = app_config.MAX_THREADS_PER_JOB

    if int(threads) == 1:
        return bam_path

    interval_files = get_region_file_list()
    for interval_index in range(len(interval_files)):
        splitted_bam_filename = FileFormats.SPLITTED_BAM_FILENAME.format(
            split_no=interval_index, name=bam_path.split(".")[0]
        )
        output_list.append(splitted_bam_filename)
        if os.path.isfile(splitted_bam_filename):
            print("BAM file already exists, skipping.")
            continue
        command = [
            "samtools",
            "view",
            "-@",
            str(app_config.MAX_THREADS_PER_JOB),
            "-b",
            bam_path,
            "-L",
            interval_files[interval_index],
            "-o",
            splitted_bam_filename,
        ]
        run(command)

        samtools_index_bam(splitted_bam_filename)

    return output_list

def samtools_index_bam(path: str):
    app_config = AppConfig()
    command = [
        "samtools",
        "index",
        "-@",
        str(app_config.MAX_THREADS_PER_JOB),
        path
    ]
    run(command)