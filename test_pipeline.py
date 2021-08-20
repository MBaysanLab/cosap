from cosap import Pipeline, FastqReader, Mapper, VariantCaller, version
from cosap import VarScanVariantCaller
from cosap import Bowtie2Mapper
from cosap import BWAMapper
from cosap._pipeline_config import PipelineKeys


germline_files = [
    FastqReader(
        "/home/mae/Desktop/sample/case_1_germline/SRR3182423_1.fastq",
        platform="illumina",
        read=1,
    ),
    FastqReader(
        "/home/mae/Desktop/sample/case_1_germline/SRR3182423_2.fastq",
        platform="illumina",
        read=2,
    ),
]

tumor_files = [
    FastqReader(
        "/home/mae/Desktop/sample/case_1_tech_rep_2_wes/SRR3182444_1.fastq",
        platform="illumina",
        read=1,
    ),
    FastqReader(
        "/home/mae/Desktop/sample/case_1_tech_rep_2_wes/SRR3182444_2.fastq",
        platform="illumina",
        read=2,
    ),
]
mapper_germline = Mapper(library="bwa", reads=germline_files, params={})
mapper_tumor = Mapper(library="bwa", reads=tumor_files, params={})
caller = VariantCaller(
    library="varscan", germline=mapper_germline, tumor=mapper_tumor, params={}
)

pipeline = Pipeline().add(mapper_germline).add(mapper_tumor).add(caller)

pipeline_config = pipeline.build()

Bowtie2Mapper().map(mapper_config=pipeline_config[PipelineKeys.MAPPING][0])

# VarScanVariantCaller().call_variants(caller_config=pipeline_config[PipelineKeys.VARIANT_CALLING][0])
