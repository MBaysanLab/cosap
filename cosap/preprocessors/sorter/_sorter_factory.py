from ._sorter import NovoalignSorter, SamtoolsSorter, _Sorter


class SorterFactory:
    SAMTOOLS_SORTER = "samtools"
    NOVOALIGN_SORTER = "novoalign"

    @classmethod
    def create(cls, processor_type: str) -> _Sorter:
        if processor_type == cls.SAMTOOLS_SORTER:
            sorter = SamtoolsSorter
        elif processor_type == cls.NOVOALIGN_SORTER:
            sorter = NovoalignSorter
        else:
            raise Exception(f"Unknown mapper type for sorter: {processor_type}.")

        return sorter
