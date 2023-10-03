import re
from abc import ABC, abstractmethod

from .._pipeline_config import *
from .._utils import convert_vcf_to_json, join_paths, read_vcf_into_df


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
            join_paths(
                self.pipeline_workdir,
                qc_dir,
                "raw_data_qualimapReport",
                "coverage_histogram.txt",
            )
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
        significant_variants = None
        uncertain_variants = None

        if "Classification" in vcf_df.columns:
            significant_variants = (
                vcf_df[
                    (vcf_df.Classification.str.contains("strong", case=False))
                    | (vcf_df.Classification.str.contains("pathogenic", case=False))
                    | (vcf_df.Classification.str.contains("potential", case=False))
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
        return join_paths(self.pipeline_workdir, vcf)


# Adapted from MultiQC https://github.com/ewels/MultiQC
def parse_qualimap_genome_results(path: str) -> dict:
    """
    Parse qualimap genome_results.txt and return metrics as dict.
    """
    file_content = open(path).read()

    results_dict = dict()
    """Parse the contents of the Qualimap BamQC genome_results.txt file"""
    regexes = {
        "bam_file": r"bam file = (.+)",
        "total_reads": r"number of reads = ([\d,]+)",
        "mapped_reads": r"number of mapped reads = ([\d,]+)",
        "mapped_bases": r"number of mapped bases = ([\d,]+)",
        "sequenced_bases": r"number of sequenced bases = ([\d,]+)",
        "mean_insert_size": r"mean insert size = ([\d,\.]+)",
        "median_insert_size": r"median insert size = ([\d,\.]+)",
        "mean_mapping_quality": r"mean mapping quality = ([\d,\.]+)",
        "general_error_rate": r"general error rate = ([\d,\.]+)",
        "mean_coverage": r"mean coverageData = ([\d,\.]+)",
    }
    d = dict()
    for k, r in regexes.items():
        r_search = re.search(r, file_content, re.MULTILINE)
        if r_search:
            if "\d" in r:
                try:
                    d[k] = float(r_search.group(1).replace(",", ""))
                except ValueError:
                    d[k] = r_search.group(1)
            else:
                d[k] = r_search.group(1)

    results_dict["total_reads"] = d["total_reads"]
    results_dict["mapped_reads"] = d["mapped_reads"]
    d["percentage_aligned"] = round((d["mapped_reads"] / d["total_reads"]) * 100, 2)
    results_dict["percentage_aligned"] = d["percentage_aligned"]
    results_dict["general_error_rate"] = d["general_error_rate"]
    results_dict["mean_coverage"] = int(d["mean_coverage"])

    return results_dict


def parse_qualimap_coverage_histogram(path) -> dict:
    """
    Parse qualimap coverage_histogram.txt and return metrics as dict.
    """

    file_content = open(path).readlines()

    results_dict = dict()
    d = dict()
    for l in file_content:
        if l.startswith("#"):
            continue
        coverage, count = l.split(None, 1)
        coverage = int(round(float(coverage)))
        count = float(count)
        d[coverage] = count

    # Find median without importing anything to do it for us
    num_counts = sum(d.values())
    cum_counts = 0
    total_cov = 0
    median_coverage = None
    for thiscov, thiscount in d.items():
        cum_counts += thiscount
        total_cov += thiscov * thiscount
        if cum_counts >= num_counts / 2:
            median_coverage = thiscov
            break

    results_dict["median_coverage"] = median_coverage

    results_dict["hist"] = d
    return results_dict
