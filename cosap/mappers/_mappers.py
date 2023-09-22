from abc import ABC, abstractmethod

from .._config import AppConfig
from .._pipeline_config import MappingKeys
import re

class _Mapper(ABC):
    @classmethod
    def _samtools_sort_command(cls, app_config: AppConfig, output_path: str):
        command = [
            "samtools",
            "sort",
            "-@",
            "6",
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

        illumina_match = re.match(r'^@([^\s]+)\s+(.+)$', header)

        readgroup = {}
        if illumina_match:
            info = illumina_match.group(1)
            instrument_id, run_number, flowcell_id, lane, tile, x, y, umi = info.split(":")

            readgroup[MappingKeys.RG_ID] = f"{instrument_id}.{flowcell_id}.{lane}"
            readgroup[MappingKeys.RG_PU] = f"{instrument_id}.{flowcell_id}.{lane}.{tile}"
            readgroup[MappingKeys.RG_PL] = "ILLUMINA"
        
        else:
            raise ValueError(f"Could not parse read group from FASTQ header: {header}")

        return readgroup
    
    @abstractmethod
    def map(self):
        pass


class _Mappable:
    pass
