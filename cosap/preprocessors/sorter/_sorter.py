import glob
import os
from abc import ABC, abstractmethod
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from .._utils import join_paths


class _Sorter(ABC):
    @abstractmethod
    def sort(self):
        pass


class _Sortable:
    @classmethod
    def _list_bam_files(cls, pipeline_config: PipelineConfig) -> List:
        return [join_paths(pipeline_config.BAM_DIR, filename) for filename in join_paths(pipeline_config.BAM_DIR, "*.bam")]


class SamtoolsSorter(_Sorter, _Sortable):
    @classmethod
    def _create_command(cls, bam_filename: str, pipeline_config: PipelineConfig, library_paths: LibraryPaths) -> str:
        command = [
            "samtools",
            "view",
            "-@",
            pipeline_config.MAPPER_THREADS,
            "-bS",
            bam_filename,
            "|",
            "samtools",
            "sort",
            "-@",
            pipeline_config.MAPPER_THREADS,,
            "-o",
            f"sorted_{bam_filename}",
        ]
        return command

    @classmethod
    def sort(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()
        bam_files = cls._list_bam_files(pipeline_config)
        for bam_filename in bam_files:
            command = cls._create_command(bam_filename=bam_filename, pipeline_config=pipeline_config, library_paths=library_paths)
            run(command, cwd=pipeline_config.BAM_DIR)
            # todo: finish this indexer
            index_command = cls._get_index_command()


class NovoalignSorter(_Sorter, _Sortable):
    @classmethod
    def _create_command(cls, bam_filename: str, pipeline_config: PipelineConfig, library_paths: LibraryPaths) -> str:
        command = [
            join_paths(library_paths.NOVOALIGN, "novosort"),
            "-m",
            "16g",
            "-t",
            ".",
            "-c",
            pipeline_config.MAPPER_THREADS,
            "--removeduplicates",
            "--keeptags",
            bam_filename,
            "-i",
            "-o",
            f"sorted_{bam_filename}",
        ]
        return command

    @classmethod
    def sort(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()
        bam_files = cls._list_bam_files(pipeline_config)
        for bam_filename in bam_files:
            command = cls._create_command(bam_filename=bam_filename, pipeline_config=pipeline_config, library_paths=library_paths)
            run(command, cwd=pipeline_config.BAM_DIR)

