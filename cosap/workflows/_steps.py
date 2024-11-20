from .._pipeline_config import PipelineKeys
from ..pipeline_builder import (MDUP, Annotator, BamReader, FastqReader,
                                Mapper, Recalibrator, Trimmer, VariantCaller)
from ..tools.annotators import AnnotatorFactory
from ..tools.mappers import MapperFactory
from ..tools.preprocessors import PreprocessorFactory
from ..tools.variant_callers import VariantCallerFactory


def trim(
    fastqs: list,
    output: str,
):

    trimmer = PreprocessorFactory().get_preprocessor("trimmer")

    readers = [FastqReader(fastq, read=i) for i, fastq in enumerate(fastqs)]
    trimmer_builder = Trimmer(
        input_step=readers,
    )
