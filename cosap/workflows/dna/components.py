from dataclasses import dataclass

from ...pipeline_builder import (
    MDUP,
    FastqReader,
    Mapper,
    Pipeline,
    Recalibrator,
    Sorter,
    Trimmer,
)
from ..models import PreprocessingOptions, SampleInput


@dataclass(frozen=True)
class ProcessedSample:
    sample: SampleInput
    reads: list
    trimmer: Trimmer
    mapper: Mapper
    mdup: MDUP
    final_bam: object


def create_fastq_readers(sample: SampleInput) -> list:
    return [
        FastqReader(sample.fastqs[0], name=sample.name, read=1),
        FastqReader(sample.fastqs[1], name=sample.name, read=2),
    ]


def preprocess_sample(
    pipeline: Pipeline,
    sample: SampleInput,
    options: PreprocessingOptions,
) -> ProcessedSample:
    reads = create_fastq_readers(sample)
    trimmer = Trimmer(input_step=reads, name=sample.name)
    pipeline.add(trimmer)

    mapper = Mapper(
        library=options.mapper,
        input_step=trimmer,
        name=options.mapper_name,
        params={
            "read_groups": {
                "ID": "0",
                "SM": sample.name,
                "PU": "0",
                "PL": options.read_group_platform,
                "LB": "0",
            }
        },
    )
    pipeline.add(mapper)

    mdup_input = mapper
    if options.sort_before_mark_duplicates:
        sorter = Sorter(input_step=mapper, name=sample.name)
        pipeline.add(sorter)
        mdup_input = sorter

    mdup = MDUP(
        input_step=mdup_input,
        name=sample.name if options.sort_before_mark_duplicates else None,
        spark=options.mark_duplicates_spark,
    )
    pipeline.add(mdup)

    final_bam = mdup
    if not options.skip_recalibration:
        recalibrator = Recalibrator(
            input_step=mdup,
            name=mdup.name,
            bed_file=options.bed_file,
        )
        pipeline.add(recalibrator)
        final_bam = recalibrator

    return ProcessedSample(
        sample=sample,
        reads=reads,
        trimmer=trimmer,
        mapper=mapper,
        mdup=mdup,
        final_bam=final_bam,
    )
