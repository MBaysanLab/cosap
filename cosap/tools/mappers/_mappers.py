import re
from abc import ABC, abstractmethod

from ..._config import AppConfig
from ..._pipeline_config import MappingKeys


class _Mapper(ABC):
    @classmethod
    def _samtools_sort_command(cls, app_config: AppConfig, output_path: str):
        command = [
            "samtools",
            "sort",
            "-@",
            "6",  # TODO: Calculate using available memory (http://www.htslib.org/doc/samtools-sort.html)
            "-n",
            "-o",
            output_path,
            "-",
        ]
        return command

    @classmethod
    def _samtools_index_command(cls, app_config: AppConfig, input_path: str):
        command = [
            "samtools",
            "index",
            input_path,
            "-@",
            "6",
        ]
        return command

    @classmethod
    def _create_readgroup_from_fastq_header(cls, fastq_file: str):
        """
        Create read groups from a FASTQ file.

        Args:
            fastq_file (str): The FASTQ file to parse.

        Returns:
            dict: A dictionary containing read group information.
        """

        def _read_fastq_header(fastq_file: str):
            if fastq_file.endswith(".gz"):
                import gzip

                with gzip.open(fastq_file, "rt") as f:
                    header = f.readline()
            else:
                with open(fastq_file, "r") as f:
                    header = f.readline()

            return header

        header = _read_fastq_header(fastq_file)

        illumina_match = re.match(
            r"^@\w+:\d+:\w+:\d+:\d+:\d+:\d+:[a-zA-Z+]+\s\d+:\w+:\d+:\w+", header
        )

        readgroup = {}
        if illumina_match:
            info = illumina_match.group(1)
            instrument_id, run_number, flowcell_id, lane, tile, x, y, umi = info.split(
                ":"
            )

            readgroup[MappingKeys.RG_ID] = f"{instrument_id}.{flowcell_id}.{lane}"
            readgroup[MappingKeys.RG_PU] = (
                f"{instrument_id}.{flowcell_id}.{lane}.{tile}"
            )
            readgroup[MappingKeys.RG_PL] = "ILLUMINA"

        return readgroup

    @classmethod
    def _create_readgroup_flags(cls, mapper_config: dict):
        """
        Take both user defined and automatically generated read group information and create a list of flags.

        Args:
            fastq_file (str): The FASTQ file to parse.

        Returns:
            dict: A dictionary containing read group flags.
        """

        automaticaly_generated_read_group = cls._create_readgroup_from_fastq_header(
            fastq_file=mapper_config[MappingKeys.INPUT]["1"]
        )

        user_defined_read_group = mapper_config[MappingKeys.PARAMS][
            MappingKeys.READ_GROUP
        ]
        return {**automaticaly_generated_read_group, **user_defined_read_group}

    @abstractmethod
    def map(self):
        pass


class _Mappable:
    pass
