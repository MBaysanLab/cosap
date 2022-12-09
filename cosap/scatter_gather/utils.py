import os
from subprocess import run
from typing import List

from .._config import AppConfig
from .._formats import FileFormats, OutputFolders
from .._library_paths import LibraryPaths
from .._utils import join_paths


def merge_bam_files(bam_list: List, output_path: str):
    command = ["samtools", "merge", output_path, *bam_list]
    run(command)


def merge_vcf_files(vcf_list: List, output_path: str):
    command = ["picard", "GatherVcfs", "-I", *vcf_list, "-O", output_path]
    run(command)


def convert_intervalfile_to_bed(interval_file_path: str):
    output_filename: interval_file_path.split(".")[0]
    command = [
        "picard",
        "IntervalListToBed",
        "-I",
        interval_file_path,
        "-O",
        f"{output_filename}.bed",
    ]
    run(command)
    os.remove(interval_file_path)


def create_gatk_intervals_as_bed(
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
        convert_intervalfile_to_bed(interval_file)


def split_bam_by_intervals(bam_path: str) -> list:
    """Splits bam by the number of threads and returns list of names of splitted bams"""

    app_config = AppConfig()
    library_paths: LibraryPaths()
    output_list = []

    threads = app_config.MAX_THREADS_PER_JOB

    if int(threads) == 1:
        return bam_path

    intervals_dir = join_paths(library_paths.INTERVALS, f"{threads}_intervals")

    if not os.path.exists(intervals_dir):
        create_gatk_intervals_as_bed(
            library_path=library_paths,
            scatter_count=threads,
            output_dir=intervals_dir,
        )

    interval_files = os.listdir(intervals_dir)
    for interval_index in range(len(interval_files)):
        splitted_bam_filename = FileFormats.SPLITTED_BAM_FILENAME.format(
            split_no=interval_index, name=bam_path.split(".")[0]
        )

        command = [
            "samtools",
            "view",
            "-b",
            bam_path,
            "-L",
            interval_files[interval_index],
            "-o",
            splitted_bam_filename,
        ]
        run(command, cwd=OutputFolders.TEMP)
        output_list.append(splitted_bam_filename)

    return output_list
