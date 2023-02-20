from itertools import groupby, product
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
        bam_qc=None,
        annotation=None,
    ):
        self.analysis_type = analysis_type
        self.normal_sample = normal_sample
        self.tumor_samples = tumor_samples
        self.normal_sample_name = normal_sample_name
        self.tumor_sample_name = tumor_sample_name
        self.bed_file = bed_file
        self.mappers = mappers
        self.variant_callers = variant_callers
        self.bam_qc = bam_qc
        self.annotation = annotation
        self.workdir = workdir
        self.pipeline = Pipeline()
        self.config = None
        self.pipeline_groups = self._get_all_pipeline_combs()

        # Do not force user for list of tuples, convert instead.
        if isinstance(self.tumor_samples, tuple):
            self.tumor_samples = [self.tumor_samples]

        if self.analysis_type.lower() == "somatic":
            self._create_somatic_pipeline()
        elif self.analysis_type.lower() == "germline":
            self._create_germline_pipeline()
        else:
            raise Exception("analysis_type can be either 'somatic' or 'germline'.")

        if "vardict" in list(map(str.lower, self.variant_callers)):
            if self.bed_file is None:
                raise Exception(
                    "Bed file should be provided for VarDict variant caller"
                )

        self._build_config()

    def _get_all_pipeline_combs(self) -> dict:
        possible_pipelines = list(product(self.mappers, self.variant_callers))
        gb_iter = groupby(possible_pipelines, lambda x: x[0])

        pipeline_groups = {}
        for key, group in gb_iter:
            key_and_group = {key: list(map(lambda x: x[1], list(group)))}
            pipeline_groups |= key_and_group

        return pipeline_groups

    def _create_somatic_pipeline(self):

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

        for i, tumor_sample in enumerate(self.tumor_samples):
            
            tumor_sample_name = self.tumor_sample_name[i]
            tumor_sample_reader = [
                FastqReader(tumor_sample[0], name=tumor_sample_name, read=1),
                FastqReader(tumor_sample[1], name=tumor_sample_name, read=2),
            ]
            trimmer_tumor = Trimmer(
                input_step=tumor_sample_reader, name=tumor_sample_name
            )
            self.pipeline.add(trimmer_tumor)

            for mapper, variant_callers in self.pipeline_groups.items():

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
                            "SM": tumor_sample_name,
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

                if self.bam_qc is not None:
                    quality_controller_tumor = QualityController(
                        library=self.bam_qc,
                        input_step=bqsr_tumor,
                        bed_file=self.bed_file,
                    )
                    self.pipeline.add(quality_controller_tumor)

                for variant_caller in variant_callers:
                    variant_caller = VariantCaller(
                        library=variant_caller,
                        germline=bqsr_normal if self.normal_sample else None,
                        tumor=bqsr_tumor,
                        bed_file=self.bed_file,
                        params={
                            "germline_sample_name": self.normal_sample_name,
                            "tumor_sample_name": tumor_sample_name,
                        },
                    )
                    self.pipeline.add(variant_caller)

                    if self.annotation is not None:
                        annotator = Annotator(
                            input_step=variant_caller, library=self.annotation
                        )
                        self.pipeline.add(annotator)

    def _create_germline_pipeline(self):

        normal_sample_reader = [
            FastqReader(self.normal_sample[0], name=self.normal_sample_name, read=1),
            FastqReader(self.normal_sample[1], name=self.normal_sample_name, read=2),
        ]
        trimmer_normal = Trimmer(
            input_step=normal_sample_reader, name=self.normal_sample_name
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

            if self.bam_qc is not None:
                quality_controller_tumor = QualityController(
                    library=self.bam_qc, input_step=bqsr_normal
                )
                self.pipeline.add(quality_controller_tumor)

            for variant_caller in variant_callers:
                variant_caller = VariantCaller(
                    library=variant_caller,
                    germline=bqsr_normal if self.normal_sample else None,
                    tumor=None,
                    bed_file=self.bed_file,
                    params={
                        "germline_sample_name": self.normal_sample_name,
                        "tumor_sample_name": self.tumor_sample_name,
                    },
                )
                self.pipeline.add(variant_caller)

                if self.annotation is not None:
                    annotator = Annotator(
                        input_step=variant_caller, library=self.annotation
                    )
                    self.pipeline.add(annotator)

    def _build_config(self):
        self.config = self.pipeline.build(workdir=self.workdir)
        return self.config

    def run_pipeline(self):
        pipeline_runner = PipelineRunner()
        pipeline_runner.run_pipeline(self.config)
        return self.config
