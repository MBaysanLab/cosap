from ..pipeline_builder import (
    FastqReader,
    Trimmer,
)
from ..tools.preprocessors import PreprocessorFactory


def trim(
    fastqs: list,
    output: str,
):
    trimmer = PreprocessorFactory().get_preprocessor("trimmer")

    readers = [FastqReader(fastq, read=i) for i, fastq in enumerate(fastqs)]
    trimmer_builder = Trimmer(
        input_step=readers,
    )
