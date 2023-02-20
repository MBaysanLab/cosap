import os
import re
from subprocess import run
from typing import List

import pandas as pd
import json


def join_paths(path: str, *paths) -> str:
    """Joins and normalizes paths according to os standards"""
    return os.path.normpath(os.path.join(path, *paths))


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


def read_vcf_into_df(path: str) -> pd.DataFrame:
    import io

    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
        df = pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#Chr": "Chr", "Ref.Gene": "Gene", "Func.refGene": "Function"})

    df.reset_index(inplace=True)
    df.rename(columns={"index": "id"}, inplace=True)
    return df


def convert_vcf_to_json(path: str) -> list:
    """
    Returns list of variants as json objects.
    """
    vcf_df = read_vcf_into_df(path)
    return vcf_df.to_dict('records')
