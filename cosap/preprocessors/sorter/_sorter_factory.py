from ._sorter import SamtoolsSorter, _Sorter


class SorterFactory:
    SAMTOOLS_SORTER = "samtools"

    @classmethod
    def create(cls, processor_type: str) -> _Sorter:
        if processor_type == cls.SAMTOOLS_SORTER:
            sorter = SamtoolsSorter
        else:
            raise Exception(f"Unknown mapper type for sorter: {processor_type}.")

        return sorter
