from ._sorter import NovoalignSorter, SamtoolsSorter, _Sorter


class SorterFactory:
    BWA_SORTER = "bwa"
    BOWTIE2_SORTER = "bowtie2"
    NOVOALIGN_SORTER = "novoalign"

    @classmethod
    def create(cls, mapper_type: str) -> _Sorter:
        mapper_type = str(mapper_type).lower()

        if mapper_type in (cls.BWA_SORTER, cls.BOWTIE2_SORTER):
            sorter = SamtoolsSorter
        elif mapper_type == cls.NOVOALIGN_SORTER:
            sorter = NovoalignSorter
        else:
            raise Exception(f"Unknown mapper type for sorter: {mapper_type}.")

        return sorter
