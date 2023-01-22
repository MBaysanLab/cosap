from abc import ABC, abstractmethod
from .._utils import (
    parse_qualimap_coverage_histogram,
    parse_qualimap_genome_results,
    convert_vcf_to_json,
    join_paths,
    read_vcf_into_df,
)
from .._pipeline_config import *


class _Parser(ABC):
    @abstractmethod
    def parse(self):
        pass


class _Parsable:
    pass


class ProjectResultsParser:
    def __init__(self, pipeline_config):
        self.pipeline_config = pipeline_config
        self.pipeline_workdir = pipeline_config[PipelineKeys.WORKDIR]
        self.vcf = self._get_vcf()
        self.qc_coverage_histogram = self._parse_qc_coverage_histogram()
        self.qc_genome_results = self._parse_qc_genome_results()
        self.variants = self._parse_vcf()
        self.variant_stats = self._parse_variant_stats()

    def _parse_qc_coverage_histogram(self):

        # If there are multiple combinations, parse the first one
        # TODO: Taking the first one might not be the preferred way.
        qc_dir = self.pipeline_config[PipelineKeys.QUALITY_CONTROL][
            list(self.pipeline_config[PipelineKeys.QUALITY_CONTROL].keys())[0]
        ][QualityControlKeys.OUTPUT]
        return parse_qualimap_coverage_histogram(
            join_paths(self.pipeline_workdir, qc_dir, "raw_data_qualimapReport", "coverage_histogram.txt")
        )

    def _parse_qc_genome_results(self):
        qc_dir = self.pipeline_config[PipelineKeys.QUALITY_CONTROL][
                list(self.pipeline_config[PipelineKeys.QUALITY_CONTROL].keys())[0]
            ][QualityControlKeys.OUTPUT]
        return parse_qualimap_genome_results(
            join_paths(self.pipeline_workdir, qc_dir, "genome_results.txt")
        )

    def _parse_vcf(self):
        return convert_vcf_to_json(join_paths(self.pipeline_workdir, self.vcf))

    def _parse_variant_stats(self):
        vcf_df = read_vcf_into_df(self.vcf)
        total_variants = vcf_df.shape[0]

        if "Classification" in vcf_df.columns:
            significant_variants = (
                vcf_df[
                    vcf_df.Classification.str.contains("strong", case=False)
                    or vcf_df.Classification.str.contains("pathogenic", case=False)
                ]
            ).shape[0]
            uncertain_variants = (
                vcf_df[vcf_df.Classification.str.contains("uncertain", case=False)]
            ).shape[0]

        return {
            "total_variants": total_variants,
            "significant_variants": significant_variants,
            "uncertain_variants": uncertain_variants,
        }

    def _get_vcf(self):
        if PipelineKeys.ANNOTATION in self.pipeline_config.keys():
            vcf = self.pipeline_config[PipelineKeys.ANNOTATION][
                list(self.pipeline_config[PipelineKeys.ANNOTATION].keys())[0]
            ][AnnotatorKeys.OUTPUT]
        elif PipelineKeys.VARIANT_CALLING in self.pipeline_config.keys():
            vcf = self.pipeline_config[PipelineKeys.VARIANT_CALLING][
                list(self.pipeline_config[PipelineKeys.VARIANT_CALLING].keys())[0]
            ][VariantCallingKeys.SNP_OUTPUT]
        return vcf
