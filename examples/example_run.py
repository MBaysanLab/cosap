from cosap import (MDUP, BamReader, DNAPipeline, DNAPipelineInput, FastqReader,
                   GeneFusionCaller, Mapper, MSICaller, Pipeline,
                   PipelineRunner, Recalibrator, Sorter, Trimmer,
                   VariantCaller)

# Reading the fastq files
normal_fastqs = [
    FastqReader("./data/SRR7890851_1.fastq.gz", name="normal", read=1),
    FastqReader("./data/SRR7890851_2.fastq.gz", name="normal", read=2),
]

tumor_fastqs = [
    FastqReader("./data/SRR7890850_1.fastq.gz", name="tumor", read=1),
    FastqReader("./data/SRR7890850_2.fastq.gz", name="tumor", read=2),
]

# Creating trimmers
trimmer_normal = Trimmer(input_step=normal_fastqs, name="normal")
trimmer_tumor = Trimmer(input_step=tumor_fastqs, name="tumor")

# Creating mappers
mapper_normal = Mapper(
    input_step=trimmer_normal,
    library="bwa",
    params={
        "read_groups": {
            "ID": "0",
            "SM": "normal",
            "PU": "0",
            "PL": "il",
            "LB": "0",
        }
    },
)
mapper_tumor = Mapper(
    input_step=trimmer_tumor,
    library="bwa",
    params={
        "read_groups": {
            "ID": "0",
            "SM": "tumor",
            "PU": "0",
            "PL": "il",
            "LB": "0",
        }
    },
)

# Creating preprocessors
mdup_normal = MDUP(input_step=mapper_normal)
mdup_tumor = MDUP(input_step=mapper_tumor)

basecal_normal = Recalibrator(input_step=mdup_normal)
basecal_tumor = Recalibrator(input_step=mdup_tumor)

# Creating variant callers
mutect = VariantCaller(library="mutect", germline=basecal_normal, tumor=basecal_tumor)

# Creating gene fusion caller and microsatellite instability caller
genefusioncaller = GeneFusionCaller(input_step=tumor_fastqs, library="genefuse")
msicaller = MSICaller(normal=basecal_normal, tumor=basecal_tumor, library="msisensor")

# Creating pipeline and adding steps to it
pipeline = (
    Pipeline()
    .add(trimmer_normal)
    .add(trimmer_tumor)
    .add(mapper_normal)
    .add(mapper_tumor)
    .add(mdup_normal)
    .add(mdup_tumor)
    .add(basecal_normal)
    .add(basecal_tumor)
    .add(mutect)
    .add(genefusioncaller)
    .add(msicaller)
)

# Creating the config that contains all information about the pipeline
config = pipeline.build()

# Create a pipeline runner and run the config file
pipeline_runner = PipelineRunner()
pipeline_runner.run_pipeline(config)


# Optionally you can just call DNAPipeline which runs all steps above on given input and algorithms.

# dna_pipeline_input = DNAPipelineInput(
#     ANALYSIS_TYPE="germline",
#     WORKDIR="/workdir/",
#     NORMAL_SAMPLE=("/data/SRR7890874_1.fastq.gz","/data/SRR7890874_2.fastq.gz"),
#     MAPPERS=["bwa2"],
#     VARIANT_CALLERS=["deepvariant"],
#     BAM_QC = "qualimap",
#     GVCF = True
# )


# dna_pipeline = DNAPipeline(dna_pipeline_input
# )
# config = dna_pipeline.run_pipeline()
