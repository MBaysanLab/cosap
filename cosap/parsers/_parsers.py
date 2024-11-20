import os
import re
from abc import ABC, abstractmethod
from glob import glob

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
        self.variants = self._get_variants()
        self.qc_coverage_histogram = self._parse_qc_coverage_histogram()
        self.qc_genome_results = self._parse_qc_genome_results()
        self.msi_score = self._parse_msi_score()

    def _parse_qc_coverage_histogram(self):
        # If there are multiple combinations, parse the first one
        # TODO: Taking the first one might not be the preferred way.

        try:
            qc_dir = self.pipeline_config[PipelineKeys.QUALITY_CONTROL][
                list(self.pipeline_config[PipelineKeys.QUALITY_CONTROL].keys())[0]
            ][QualityControlKeys.OUTPUT]
        except IndexError:
            qc_dir = None

        if qc_dir is None:
            return None

        return parse_qualimap_coverage_histogram(
            join_paths(
                self.pipeline_workdir,
                qc_dir,
                "raw_data_qualimapReport",
                "coverage_histogram.txt",
            )
        )

    def _parse_qc_genome_results(self):
        try:
            qc_dir = self.pipeline_config[PipelineKeys.QUALITY_CONTROL][
                list(self.pipeline_config[PipelineKeys.QUALITY_CONTROL].keys())[0]
            ][QualityControlKeys.OUTPUT]
        except (KeyError, IndexError):
            qc_dir = None

        if qc_dir is None:
            return None

        return parse_qualimap_genome_results(
            join_paths(self.pipeline_workdir, qc_dir, "genome_results.txt")
        )

    def _parse_vcf(self, vcf, caller_type, sample_name):
        return convert_vcf_to_json(
            join_paths(self.pipeline_workdir, vcf), caller_type, sample_name
        )

    def _get_variants(self):

        # If there is variant caller in the pipeline config, parse the VCF
        if PipelineKeys.VARIANT_CALLING in self.pipeline_config:
            variant_caller_config = self.pipeline_config[PipelineKeys.VARIANT_CALLING][
                list(self.pipeline_config[PipelineKeys.VARIANT_CALLING].keys())[0]
            ]
            variant_caller = variant_caller_config[VariantCallingKeys.LIBRARY]
            vcf = variant_caller_config[VariantCallingKeys.ALL_VARIANTS_OUTPUT]
        else:
            raise ValueError("No variant caller found in pipeline config.")

        # If tumor sample is present, parse it, otherwise parse the normal sample
        if VariantCallingKeys.TUMOR_INPUT in variant_caller_config:
            sample_name = variant_caller_config[VariantCallingKeys.PARAMS][
                VariantCallingKeys.TUMOR_SAMPLE_NAME
            ]
        elif VariantCallingKeys.GERMLINE_INPUT in variant_caller_config:
            sample_name = variant_caller_config[VariantCallingKeys.PARAMS][
                VariantCallingKeys.GERMLINE_SAMPLE_NAME
            ]
        else:
            raise ValueError("No VCF found in pipeline config.")

        vcf_dir = variant_caller_config[VariantCallingKeys.OUTPUT_DIR]
        vcf_path = join_paths(self.pipeline_workdir, vcf_dir, vcf)

        return self._parse_vcf(
            vcf_path, caller_type=variant_caller, sample_name=sample_name
        )

    def _parse_msi_score(self):
        try:
            msi_config = self.pipeline_config[PipelineKeys.MSI][
                list(self.pipeline_config[PipelineKeys.MSI].keys())[0]
            ]
            msi_file = msi_config[MSICallingKeys.OUTPUT]
        except (IndexError, KeyError):
            return None

        return parse_msi_results(join_paths(self.pipeline_workdir, msi_file))


# Adapted from MultiQC https://github.com/ewels/MultiQC
def parse_qualimap_genome_results(path: str) -> dict:
    """
    Parse qualimap genome_results.txt and return metrics as dict.
    """

    if not os.path.exists(path):
        return None

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

    if not os.path.exists(path):
        return None

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


def parse_msi_results(path):
    """
    Read msisensor results file and return msi score.
    """

    if not os.path.exists(path):
        return None

    # The results are in the second line of the file
    with open(path) as f:
        next(f)
        line = next(f)

        msi_score = line.split()[2]

    return msi_score
