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
    def _get_file_information(cls, fastq_list: List, sample_type: str, trimmed: bool = False) -> Dict:
        splitted_files = list(map(lambda x: x.replace(".fastq.gz", "").split("_"), fastq_list))

        # Picks the correct dict keys in order
        if sample_type == "tumor":
            keys = ("Sample_ID", "Index", "Lanes", "Pairs", "Number_of_seq")
        elif sample_type in ("germline", "normal"):
            keys = ("Sample_ID", "Germline", "Index", "Lanes", "Pairs", "Number_of_seq")
        else:
            raise Exception(f"Unknown sample type: {sample_type}")

        # converts True to 1 and visa versa
        i = int(trimmed)
        
        # Cleans filename string into pieces
        values = [list(set(items[i:len(keys)])) for items in zip(*splitted_files)]

        # creates a dict by the order
        info_dict = dict(zip(keys, values))

        return info_dict

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




