from dataclasses import dataclass, field
from itertools import groupby, product
from typing import List, Tuple

from ..pipeline_builder import *
from ..pipeline_runner import PipelineRunner


@dataclass
class DNAPipelineInput:
    ANALYSIS_TYPE: str
    WORKDIR: str
    NORMAL_SAMPLE: tuple[str, str] = None
    TUMOR_SAMPLES: List[Tuple[str, str]] = field(default_factory=lambda: None)
    BED_FILE: str = None
    MAPPERS: List[str] = field(default_factory=lambda: ["bwa"])
    VARIANT_CALLERS: List[str] = field(default_factory=lambda: [""])
    NORMAL_SAMPLE_NAME: str = "normal"
    TUMOR_SAMPLE_NAME: str = "tumor"
    BAM_QC: str = None
    ANNOTATORS: List[str] = field(default_factory=lambda: None)
    GVCF: bool = False
    MSI: bool = False
    GENEFUSION: bool = False

    def __post_init__(self):
        if isinstance(self.TUMOR_SAMPLES, tuple):
            self.TUMOR_SAMPLES = [self.TUMOR_SAMPLES]

        if isinstance(self.ANNOTATORS, str):
            self.ANNOTATORS = [self.ANNOTATORS]

        if "vardict" in list(map(str.lower, self.VARIANT_CALLERS)):
            if self.BED_FILE is None:
                raise Exception(
                    "Bed file should be provided for VarDict variant caller"
                )


class DNAPipeline:
    def __init__(
        self,
        dna_pipeline_input: DNAPipelineInput,
    ):
        self.input = dna_pipeline_input
        self.pipeline = Pipeline()
        self.config = None
        self.pipeline_groups = self._get_all_pipeline_combs()

        if self.input.ANALYSIS_TYPE.lower() == "somatic":
            self._create_somatic_pipeline()
        elif self.input.ANALYSIS_TYPE.lower() == "germline":
            self._create_germline_pipeline()
        else:
            raise Exception(
                "input.ANALYSIS_TYPE can be either 'somatic' or 'germline'."
            )

        self._build_config()

    def _get_all_pipeline_combs(self) -> dict:

        """
        Returns all possible combinations of mappers and variant callers
        """

        possible_pipelines = list(
            product(self.input.MAPPERS, self.input.VARIANT_CALLERS)
        )
        gb_iter = groupby(possible_pipelines, lambda x: x[0])

        pipeline_groups = {}
        for key, group in gb_iter:
            key_and_group = {key: list(map(lambda x: x[1], list(group)))}
            pipeline_groups |= key_and_group

        return pipeline_groups

    def _create_somatic_pipeline(self):

        if self.input.NORMAL_SAMPLE is not None:
            normal_sample_reader = [
                FastqReader(
                    self.input.NORMAL_SAMPLE[0],
                    name=self.input.NORMAL_SAMPLE_NAME,
                    read=1,
                ),
                FastqReader(
                    self.input.NORMAL_SAMPLE[1],
                    name=self.input.NORMAL_SAMPLE_NAME,
                    read=2,
                ),
            ]
            trimmer_normal = Trimmer(
                input_step=normal_sample_reader, name=self.input.NORMAL_SAMPLE_NAME
            )
            self.pipeline.add(trimmer_normal)

        for i, tumor_sample in enumerate(self.input.TUMOR_SAMPLES):

            tumor_sample_name = (
                self.input.TUMOR_SAMPLE_NAME[i]
                if isinstance(self.input.TUMOR_SAMPLE_NAME, list)
                else self.input.TUMOR_SAMPLE_NAME
            )
            tumor_sample_reader = [
                FastqReader(tumor_sample[0], name=tumor_sample_name, read=1),
                FastqReader(tumor_sample[1], name=tumor_sample_name, read=2),
            ]

            if self.input.GENEFUSION:
                genefusioncaller = GeneFusionCaller(
                    input_step=tumor_sample_reader, library="genefuse"
                )
                self.pipeline.add(genefusioncaller)

            trimmer_tumor = Trimmer(
                input_step=tumor_sample_reader, name=tumor_sample_name
            )
            self.pipeline.add(trimmer_tumor)

            for mapper, variant_callers in self.pipeline_groups.items():

                # If normal sample suplied, add related steps to pipeline
                if self.input.NORMAL_SAMPLE:
                    mapper_normal = Mapper(
                        library=mapper,
                        input_step=trimmer_normal,
                        params={
                            "read_groups": {
                                "ID": "0",
                                "SM": self.input.NORMAL_SAMPLE_NAME,
                                "PU": "0",
                                "PL": "il",
                                "LB": "0",
                            }
                        },
                    )
                    mdup_normal = MDUP(input_step=mapper_normal)
                    bqsr_normal = Recalibrator(
                        input_step=mdup_normal, bed_file=self.input.BED_FILE
                    )

                    self.pipeline.add(mapper_normal)
                    self.pipeline.add(mdup_normal)
                    self.pipeline.add(bqsr_normal)

                mapper_tumor = Mapper(
                    library=mapper,
                    input_step=trimmer_tumor,
                    params={
                        "read_groups": {
                            "ID": "0",
                            "SM": tumor_sample_name,
                            "PU": "0",
                            "PL": "il",
                            "LB": "0",
                        }
                    },
                )
                mdup_tumor = MDUP(input_step=mapper_tumor)

                if self.input.MSI:
                    msicaller = MSICaller(
                        normal=mdup_normal, tumor=mdup_tumor, library="msisensor"
                    )
                    self.pipeline.add(msicaller)

                bqsr_tumor = Recalibrator(
                    input_step=mdup_tumor, bed_file=self.input.BED_FILE
                )

                self.pipeline.add(mapper_tumor)
                self.pipeline.add(mdup_tumor)

                self.pipeline.add(bqsr_tumor)

                if self.input.BAM_QC is not None:
                    quality_controller_tumor = QualityController(
                        library=self.input.BAM_QC,
                        input_step=bqsr_tumor,
                        bed_file=self.input.BED_FILE,
                    )
                    self.pipeline.add(quality_controller_tumor)

                for variant_caller in variant_callers:
                    variant_caller = VariantCaller(
                        library=variant_caller,
                        germline=bqsr_normal if self.input.NORMAL_SAMPLE else None,
                        tumor=bqsr_tumor,
                        bed_file=self.input.BED_FILE,
                        params={
                            "germline_sample_name": self.input.NORMAL_SAMPLE_NAME,
                            "tumor_sample_name": tumor_sample_name,
                        },
                    )
                    self.pipeline.add(variant_caller)

                    if self.input.ANNOTATORS is not None:
                        for annotator in self.input.ANNOTATORS:
                            ann = Annotator(
                                input_step=variant_caller, library=annotator
                            )
                            self.pipeline.add(ann)

    def _create_germline_pipeline(self):

        normal_sample_reader = [
            FastqReader(
                self.input.NORMAL_SAMPLE[0], name=self.input.NORMAL_SAMPLE_NAME, read=1
            ),
            FastqReader(
                self.input.NORMAL_SAMPLE[1], name=self.input.NORMAL_SAMPLE_NAME, read=2
            ),
        ]
        trimmer_normal = Trimmer(
            input_step=normal_sample_reader, name=self.input.NORMAL_SAMPLE_NAME
        )
        self.pipeline.add(trimmer_normal)

        for mapper, variant_callers in self.pipeline_groups.items():
            # If normal sample suplied, add related steps to pipeline

            mapper_normal = Mapper(
                library=mapper,
                input_step=trimmer_normal,
                params={
                    "read_groups": {
                        "ID": "0",
                        "SM": self.input.NORMAL_SAMPLE_NAME,
                        "PU": "0",
                        "PL": "il",
                        "LB": "0",
                    }
                },
            )
            mdup_normal = MDUP(input_step=mapper_normal)
            bqsr_normal = Recalibrator(
                input_step=mdup_normal, bed_file=self.input.BED_FILE
            )

            self.pipeline.add(mapper_normal)
            self.pipeline.add(mdup_normal)
            self.pipeline.add(bqsr_normal)

            if self.input.BAM_QC is not None:
                quality_controller_tumor = QualityController(
                    library=self.input.BAM_QC, input_step=bqsr_normal
                )
                self.pipeline.add(quality_controller_tumor)

            for variant_caller in variant_callers:
                variant_caller = VariantCaller(
                    library=variant_caller,
                    germline=bqsr_normal if self.input.NORMAL_SAMPLE else None,
                    tumor=None,
                    bed_file=self.input.BED_FILE,
                    gvcf=self.input.GVCF,
                )
                self.pipeline.add(variant_caller)

                if self.input.ANNOTATORS is not None:
                    for annotator in self.input.ANNOTATORS:
                        ann = Annotator(input_step=variant_caller, library=annotator)
                        self.pipeline.add(ann)

    def _build_config(self):
        self.config = self.pipeline.build(workdir=self.input.WORKDIR)
        return self.config

    def run_pipeline(self):
        pipeline_runner = PipelineRunner()
        pipeline_runner.run_pipeline(self.config)
