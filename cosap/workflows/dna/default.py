from dataclasses import dataclass, field
from typing import List, Tuple

from ...pipeline_builder import (
    Annotator,
    GeneFusionCaller,
    MitochondrialVariantCaller,
    MSICaller,
    Pipeline,
    QualityController,
    VariantCaller,
)
from ...runners.pipeline_runner import PipelineRunner
from ..models import PreprocessingOptions, SampleInput
from .components import ProcessedSample, create_fastq_readers, preprocess_sample

STRUCTURAL_VARIANT_CALLERS = ["manta"]
SV_ANNOTATORS = ["annotsv"]


@dataclass
class DNAPipelineInput:
    ANALYSIS_TYPE: str
    WORKDIR: str
    NORMAL_SAMPLE: tuple[str, str] = None
    TUMOR_SAMPLES: List[Tuple[str, str]] = field(default_factory=lambda: None)
    BED_FILE: str = None
    MAPPERS: List[str] = field(default_factory=lambda: ["bwa"])
    VARIANT_CALLERS: List[str] = field(default_factory=lambda: [""])
    NORMAL_SAMPLE_NAME: str = None
    TUMOR_SAMPLE_NAME: str = None
    BAM_QC: str = None
    ANNOTATORS: List[str] = field(default_factory=lambda: None)
    GVCF: bool = False
    MSI: bool = False
    GENEFUSION: bool = False
    MITOCHONDRIAL: bool = True
    DEVICE: str = "cpu"
    SKIP_RECALIBRATION: bool = True
    MODEL_TYPE: str = None

    def __post_init__(self):
        if isinstance(self.TUMOR_SAMPLES, tuple):
            self.TUMOR_SAMPLES = [self.TUMOR_SAMPLES]

        if isinstance(self.ANNOTATORS, str):
            self.ANNOTATORS = [self.ANNOTATORS]

        if "vardict" in list(map(str.lower, self.VARIANT_CALLERS)):
            if self.BED_FILE is None:
                raise Exception("Bed file should be provided for VarDict variant caller")

        analysis = self.ANALYSIS_TYPE.lower()
        if analysis == "germline" and self.NORMAL_SAMPLE is None:
            raise ValueError(
                "A normal sample (--normal_sample) is required for germline analysis. "
                "Provide a comma-separated pair of R1,R2 FASTQ paths."
            )
        if analysis == "somatic":
            if self.NORMAL_SAMPLE is None:
                raise ValueError(
                    "A normal sample (--normal_sample) is required for somatic analysis. "
                    "Provide a comma-separated pair of R1,R2 FASTQ paths."
                )
            if not self.TUMOR_SAMPLES:
                raise ValueError(
                    "At least one tumor sample (--tumor_sample) is required for somatic analysis. "
                    "Provide a comma-separated pair of R1,R2 FASTQ paths."
                )

        if self.TUMOR_SAMPLE_NAME is None:
            self.TUMOR_SAMPLE_NAME = "tumor"

        if self.NORMAL_SAMPLE_NAME is None:
            self.NORMAL_SAMPLE_NAME = "normal"


class DNAPipeline:
    def __init__(self, dna_pipeline_input: DNAPipelineInput):
        self.input = dna_pipeline_input
        self.pipeline = self.build()
        self.config = self.pipeline.build(workdir=self.input.WORKDIR)

    def build(self) -> Pipeline:
        pipeline = Pipeline()
        if self.input.ANALYSIS_TYPE.lower() == "somatic":
            self._create_somatic_pipeline(pipeline)
        elif self.input.ANALYSIS_TYPE.lower() == "germline":
            self._create_germline_pipeline(pipeline)
        else:
            raise Exception("input.ANALYSIS_TYPE can be either 'somatic' or 'germline'.")
        return pipeline

    @property
    def pipeline_groups(self) -> dict:
        return {mapper: list(self.input.VARIANT_CALLERS) for mapper in self.input.MAPPERS}

    def _preprocess(self, pipeline: Pipeline, sample: SampleInput, mapper: str):
        return preprocess_sample(
            pipeline,
            sample,
            PreprocessingOptions(
                mapper=mapper,
                skip_recalibration=self.input.SKIP_RECALIBRATION,
                bed_file=self.input.BED_FILE,
            ),
        )

    def _add_mitochondrial_variant_caller(
        self,
        pipeline: Pipeline,
        sample: ProcessedSample,
    ):
        if self.input.MITOCHONDRIAL:
            pipeline.add(MitochondrialVariantCaller(input_bam=sample.final_bam))

    def _add_quality_control(self, pipeline: Pipeline, sample: ProcessedSample):
        if self.input.BAM_QC is not None:
            pipeline.add(
                QualityController(
                    library=self.input.BAM_QC,
                    input_step=sample.final_bam,
                    bed_file=self.input.BED_FILE,
                )
            )

    def _add_annotations(self, pipeline: Pipeline, variant_step, variant_caller: str):
        if self.input.ANNOTATORS is None:
            return

        for annotator in self.input.ANNOTATORS:
            if annotator in SV_ANNOTATORS and variant_caller not in STRUCTURAL_VARIANT_CALLERS:
                continue
            pipeline.add(Annotator(input_step=variant_step, library=annotator))

    def _add_germline_variant_callers(
        self,
        pipeline: Pipeline,
        sample: ProcessedSample,
        variant_callers: list,
    ):
        for variant_caller in variant_callers:
            vc = VariantCaller(
                library=variant_caller,
                normal_sample=sample.final_bam,
                tumor_sample=None,
                bed_file=self.input.BED_FILE,
                gvcf=self.input.GVCF,
                model_type=self.input.MODEL_TYPE,
            )
            pipeline.add(vc)
            self._add_annotations(pipeline, vc, variant_caller)

    def _add_somatic_variant_callers(
        self,
        pipeline: Pipeline,
        normal: ProcessedSample,
        tumor: ProcessedSample,
        tumor_name: str,
        variant_callers: list,
    ):
        for variant_caller in variant_callers:
            vc = VariantCaller(
                library=variant_caller,
                normal_sample=normal.final_bam,
                tumor_sample=tumor.final_bam,
                bed_file=self.input.BED_FILE,
                model_type=self.input.MODEL_TYPE,
                params={
                    "germline_sample_name": self.input.NORMAL_SAMPLE_NAME,
                    "tumor_sample_name": tumor_name,
                },
            )
            pipeline.add(vc)
            self._add_annotations(pipeline, vc, variant_caller)

    def _tumor_sample_name(self, index: int) -> str:
        if isinstance(self.input.TUMOR_SAMPLE_NAME, list):
            return self.input.TUMOR_SAMPLE_NAME[index]
        if len(self.input.TUMOR_SAMPLES) == 1:
            return self.input.TUMOR_SAMPLE_NAME
        return f"{self.input.TUMOR_SAMPLE_NAME}_{index + 1}"

    def _create_somatic_pipeline(self, pipeline: Pipeline):
        normal_input = SampleInput(
            name=self.input.NORMAL_SAMPLE_NAME,
            fastqs=self.input.NORMAL_SAMPLE,
        )
        normal_by_mapper = {
            mapper: self._preprocess(pipeline, normal_input, mapper)
            for mapper in self.input.MAPPERS
        }
        for normal in normal_by_mapper.values():
            self._add_mitochondrial_variant_caller(pipeline, normal)

        for i, tumor_fastqs in enumerate(self.input.TUMOR_SAMPLES):
            tumor_name = self._tumor_sample_name(i)
            tumor_input = SampleInput(name=tumor_name, fastqs=tumor_fastqs)
            tumor_reads = create_fastq_readers(tumor_input)

            if self.input.GENEFUSION:
                pipeline.add(GeneFusionCaller(input_step=tumor_reads, library="genefuse"))

            for mapper, variant_callers in self.pipeline_groups.items():
                normal = normal_by_mapper[mapper]
                tumor = self._preprocess(pipeline, tumor_input, mapper)

                if self.input.MSI:
                    pipeline.add(
                        MSICaller(
                            normal=normal.mdup,
                            tumor=tumor.mdup,
                            library="msisensor",
                        )
                    )

                self._add_mitochondrial_variant_caller(pipeline, tumor)
                self._add_quality_control(pipeline, tumor)
                self._add_somatic_variant_callers(
                    pipeline,
                    normal,
                    tumor,
                    tumor_name,
                    variant_callers,
                )

    def _create_germline_pipeline(self, pipeline: Pipeline):
        normal_input = SampleInput(
            name=self.input.NORMAL_SAMPLE_NAME,
            fastqs=self.input.NORMAL_SAMPLE,
        )

        for mapper, variant_callers in self.pipeline_groups.items():
            normal = self._preprocess(pipeline, normal_input, mapper)
            self._add_mitochondrial_variant_caller(pipeline, normal)
            self._add_quality_control(pipeline, normal)
            self._add_germline_variant_callers(pipeline, normal, variant_callers)

    def _build_config(self):
        self.config = self.pipeline.build(workdir=self.input.WORKDIR)
        return self.config

    def run_pipeline(self):
        pipeline_runner = PipelineRunner(device=self.input.DEVICE)
        return pipeline_runner.run_pipeline(self.config)
