from abc import ABC, abstractmethod
from .._utils import parse_qualimap_coverage_histogram, parse_qualimap_genome_results, convert_vcf_to_json


class _Parser(ABC):
    @abstractmethod
    def parse(self):
        pass


class _Parsable:
    pass


class ProjectResultsParser:
    def __init__(self, pipeline_config):
        self.pipeline_config = pipeline_config

        self.qc_coverage_histogram = self._parse_qc_coverage_histogram()
        self.qc_genome_results = self._parse_qc_genome_results()
        self.variants = self._parse_annotated_vcf()
        self.variant_stats = self._parse_variant_stats()

    
    def _parse_qc_coverage_histogram(self):
        pass

    def _parse_qc_genome_results(self):
        pass

    def _parse_annotated_vcf(self):
        pass

    def _parse_variant_stats(self):
        pass
