from ._bowtie_mapper import Bowtie2Mapper
from ._bwa_mapper import BWAMapper
from ._mappers import _Mapper


class MapperFactory:
    BWA_MAPPER = "bwa"
    BOWTIE2_MAPPER = "bowtie2"

    @classmethod
    def create(cls, mapper_type: str) -> _Mapper:
        mapper_type = str(mapper_type).lower()

        if mapper_type == cls.BWA_MAPPER:
            mapper = BWAMapper
        elif mapper_type == cls.BOWTIE2_MAPPER:
            mapper = Bowtie2Mapper
        else:
            raise Exception(f"Unknown mapper type: {mapper_type}")

        return mapper
