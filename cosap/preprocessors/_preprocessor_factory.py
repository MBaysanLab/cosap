from ._preprocessors import _Preprocessor
from ._mark_duplicate import MarkDuplicate
from ._base_recalibrator import BaseRecalibrator
from ._indexer import BamIndexer
from ._sorter import SamtoolsSorter
from ._merger import BamMerger


class PreprocessorFactory:
    markduplicate = "mark_duplicate"
    indexer = "indexer"
    base_recalibrator = "base_recalibrator"
    merger = "merger"
    sorter = "sorter"

    @classmethod
    def create(cls, preprocessor_type: str) -> _Preprocessor:
        preprocessor_type = str(preprocessor_type).lower()

        if preprocessor_type == cls.markduplicate:
            preprocessor = MarkDuplicate
        elif preprocessor_type == cls.indexer:
            preprocessor = BamIndexer
        elif preprocessor_type == cls.base_recalibrator:
            preprocessor = BaseRecalibrator
        elif preprocessor_type == cls.merger:
            preprocessor = BamMerger
        elif preprocessor_type == cls.sorter:
            preprocessor = SamtoolsSorter
        else:
            raise Exception(f"Unknown mapper type: {preprocessor_type}")

        return preprocessor
