from itertools import product
from typing import List, Tuple

from ..pipeline_builder import *
from ..pipeline_runner import PipelineRunner


class DNAPipeline:
    def __init__(
        self,
        analysis_type: str,
        workdir: str,
        normal_sample: tuple[str, str] = None,
        tumor_samples: List[Tuple[str, str]] = None,
        bed_file: str = None,
        mappers: List[str] = ["bwa"],
        variant_callers: List[str] = [""],
        normal_sample_name="normal",
        tumor_sample_name="tumor",
    ):
        self.analysis_type = analysis_type
        self.normal_sample = normal_sample
        self.tumor_samples = tumor_samples
        self.normal_sample_name = normal_sample_name
        self.tumor_sample_name = tumor_sample_name
        self.bed_file = bed_file
        self.mappers = mappers
        self.variant_callers = variant_callers
        self.workdir = workdir
        self.pipeline = Pipeline()
        self.config = None

        # Do not force user for list of tuples, convert instead.
        if isinstance(self.tumor_samples, tuple):
            self.tumor_samples = [self.tumor_samples]

        if self.analysis_type.lower() == "somatic":
            self._create_somatic_pipeline()
        elif self.analysis_type.lower() == "germline":
            self._create_germline_pipeline()
        else:
            raise Exception("analysis_type can be either 'somatic' or 'germline'.")

        self._build_config()

    def _create_somatic_pipeline(self):
        possible_pipelines = list(product(self.mappers, self.variant_callers))

        if self.normal_sample is not None:

            normal_sample_reader = [
                FastqReader(
                    self.normal_sample[0], name=self.normal_sample_name, read=1
                ),
                FastqReader(
                    self.normal_sample[1], name=self.normal_sample_name, read=2
                ),
            ]
            trimmer_normal = Trimmer(
                input_step=normal_sample_reader, name=self.normal_sample_name
            )
            self.pipeline.add(trimmer_normal)

        for tumor_sample in self.tumor_samples:

            tumor_sample_reader = [
                FastqReader(tumor_sample[0], name=self.tumor_sample_name, read=1),
                FastqReader(tumor_sample[1], name=self.tumor_sample_name, read=2),
            ]
            trimmer_tumor = Trimmer(
                input_step=tumor_sample_reader, name=self.tumor_sample_name
            )
            self.pipeline.add(trimmer_tumor)

            for pipeline_perm in possible_pipelines:
                mapper = pipeline_perm[0]
                variant_caller = pipeline_perm[1]

                # If normal sample suplied, add related steps to pipeline
                if self.normal_sample:
                    mapper_normal = Mapper(
                        library=mapper,
                        input_step=trimmer_normal,
                        params={
                            "read_groups": {
                                "ID": "0",
                                "SM": self.normal_sample_name,
                                "PU": "0",
                                "PL": "il",
                                "LB": "0",
                            }
                        },
                    )
                    mdup_normal = MDUP(input_step=mapper_normal)
                    bqsr_normal = Recalibrator(input_step=mdup_normal)

                    self.pipeline.add(mapper_normal)
                    self.pipeline.add(mdup_normal)
                    self.pipeline.add(bqsr_normal)

                mapper_tumor = Mapper(
                    library=mapper,
                    input_step=trimmer_tumor,
                    params={
                        "read_groups": {
                            "ID": "0",
                            "SM": self.tumor_sample_name,
                            "PU": "0",
                            "PL": "il",
                            "LB": "0",
                        }
                    },
                )
                mdup_tumor = MDUP(input_step=mapper_tumor)
                bqsr_tumor = Recalibrator(input_step=mdup_tumor)

                self.pipeline.add(mapper_tumor)
                self.pipeline.add(mdup_tumor)
                self.pipeline.add(bqsr_tumor)

                variant_caller = VariantCaller(
                    library=variant_caller,
                    germline=bqsr_normal if self.normal_sample else None,
                    tumor=bqsr_tumor,
                    params={"germline_sample_name": self.normal_sample_name},
                )

                annotator = Annotator(input_step=variant_caller, library="cancervar")

                self.pipeline.add(variant_caller)
                self.pipeline.add(annotator)

                quality_controller_tumor = QualityController(
                    library="qualimap", input_step=bqsr_tumor
                )
                self.pipeline.add(quality_controller_tumor)

    def _create_germline_pipeline(self):
        possible_pipelines = list(product(self.mappers, self.variant_callers))

        normal_sample_reader = [
            FastqReader(self.normal_sample[0], name=self.normal_sample_name, read=1),
            FastqReader(self.normal_sample[1], name=self.normal_sample_name, read=2),
        ]
        trimmer_normal = Trimmer(
            input_step=normal_sample_reader, name=self.normal_sample_name
        )
        self.pipeline.add(trimmer_normal)

        for pipeline_perm in possible_pipelines:
            mapper = pipeline_perm[0]
            variant_caller = pipeline_perm[1]

            # If normal sample suplied, add related steps to pipeline

            mapper_normal = Mapper(
                library=mapper,
                input_step=trimmer_normal,
                params={
                    "read_groups": {
                        "ID": "0",
                        "SM": self.normal_sample_name,
                        "PU": "0",
                        "PL": "il",
                        "LB": "0",
                    }
                },
            )
            mdup_normal = MDUP(input_step=mapper_normal)
            bqsr_normal = Recalibrator(input_step=mdup_normal)

            self.pipeline.add(mapper_normal)
            self.pipeline.add(mdup_normal)
            self.pipeline.add(bqsr_normal)

            variant_caller = VariantCaller(
                library=variant_caller,
                germline=bqsr_normal if self.normal_sample else None,
                tumor=None,
                params={"germline_sample_name": self.normal_sample_name},
            )

            annotator = Annotator(input_step=variant_caller, library="intervar")

            self.pipeline.add(variant_caller)
            self.pipeline.add(annotator)

            quality_controller = QualityController(input_step=quality_controller)
            self.pipeline.add(quality_controller)

    def _build_config(self):
        self.config = self.pipeline.build(workdir=self.workdir)
        return self.config

    def run_pipeline(self):
        pipeline_runner = PipelineRunner()
        pipeline_runner.run_pipeline(self.config)
