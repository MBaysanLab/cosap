from abc import ABC, abstractmethod
from typing import Dict, List
import glob
from .._library_paths import LibraryPaths


class _Mapper(ABC):
    @abstractmethod
    def map(self):
        pass


class _Mappable:
    @classmethod
    def _list_fastq_files(cls, file_path: str) -> List:
        return glob.glob(os.path.join(file_path, "*fastq.gz"))

    @classmethod
    def _get_file_information(
        cls, fastq_list: List, sample_type: str, trimmed: bool = False
    ) -> Dict:
        format_string = ""
        if trimmed:
            format_string = "{Drop}_"

        if sample_type == "tumor":
            format_string += "{Sample_ID}_{Index}_{Lanes}_{Pairs_Number_of_seq}.fastq.gz"
        elif sample_type in ("germline", "normal"):
            format_string += "{Sample_ID}_{Germline}_{Index}_{Lanes}_{Pairs_Number_of_seq}.fastq.gz"
        else:
            raise Exception(f"Unknown sample type: {sample_type}")

        fastq_filename_parser = Parser(format_string)

        info_dict = defaultdict(set)
        for filename in fastq_list:
            file_info = fastq_filename_parser.parse(filename)
            for key, value in file_info.items():
                if key != "Drop":
                    info_dict[key].add(value)

        return dict(map(lambda x: (x[0], list(x[1]), info_dict.items())))



class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def map(cls, pipeline_config: Dict):
        library_paths = LibraryPaths()
        pipeline_config


class Bowtie2Mapper(_Mapper, _Mappable):
    @classmethod
    def map(cls, pipeline_config: Dict):
        pass


class NovoalignMapper(_Mapper, _Mappable):
    @classmethod
    def map(cls, pipeline_config: Dict):
        pass
