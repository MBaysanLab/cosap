from ._bowtie_mapper import Bowtie2Mapper
from ._bwa2_mapper import BWA2Mapper
from ._bwa_mapper import BWAMapper
from ._mappers import _Mapper


class MapperFactory:
    BWA_MAPPER = "bwa"
    BWA2_MAPPER = "bwa2"
    BOWTIE2_MAPPER = ("bowtie", "bowtie2")

    @classmethod
    def create(cls, mapper_type: str) -> _Mapper:
        mapper_type = str(mapper_type).lower()

        if mapper_type == cls.BWA_MAPPER:
            mapper = BWAMapper
        elif mapper_type in cls.BOWTIE2_MAPPER:
            mapper = Bowtie2Mapper
        elif mapper_type == cls.BWA2_MAPPER:
            mapper = BWA2Mapper
        else:
            raise Exception(f"Unknown mapper type: {mapper_type}")

        return mapper
