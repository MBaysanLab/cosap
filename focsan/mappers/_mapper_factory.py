from ._bowtie_mapper import Bowtie2Mapper
from ._bwa_mapper import BWAMapper
from ._mappers import _Mapper
from ._novoalign_mapper import NovoalignMapper


class MapperFactory:
    BWA_MAPPER = "bwa"
    BOWTIE2_MAPPER = "bowtie2"
    NOVOALIGN_MAPPER = "novoalign"

    @classmethod
    def create(cls, mapper_type: str) -> _Mapper:
        mapper_type = str(mapper_type).lower()

        if mapper_type == cls.BWA_MAPPER:
            mapper = BWAMapper
        elif mapper_type == cls.BOWTIE2_MAPPER:
            mapper = Bowtie2Mapper
        elif mapper_type == cls.NOVOALIGN_MAPPER:
            mapper = NovoalignMapper
        else:
            raise Exception(f"Unknown mapper type: {mapper_type}")

        return mapper
