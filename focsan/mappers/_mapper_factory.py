from ._mappers import _Mapper, BWAMapper, Bowtie2Mapper


class MapperFactory:
    BWA_MAPPER = "bwa"
    BOWTIE2_MAPPER = "bowtie2"

    @classmethod
    def create(cls, mapper_type: str) -> _Mapper:
        mapper_type = mapper_type.lower()

        if mapper_type == cls.BWA_MAPPER:
            mapper = BWAMapper
        else mapper_type == cls.BOWTIE2_MAPPER:
            mapper = Bowtie2Mapper
        else:
            raise Exception(f"Unknown mapper type: {mapper_type}")

        return mapper
