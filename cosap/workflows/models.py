from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class SampleInput:
    name: str
    fastqs: Tuple[str, str]


@dataclass(frozen=True)
class PreprocessingOptions:
    mapper: str = "bwa"
    skip_recalibration: bool = True
    bed_file: Optional[str] = None
    read_group_platform: str = "ILLUMINA"
    sort_before_mark_duplicates: bool = False
    mark_duplicates_spark: bool = True
    mapper_name: Optional[str] = None


@dataclass(frozen=True)
class MitochondrialOptions:
    enabled: bool = True
    contig_name: str = "chrM"
    median_autosomal_coverage: Optional[float] = None
    min_allele_fraction: Optional[float] = None
    f_score_beta: Optional[float] = None
    contamination_estimate: Optional[float] = None
    run_haplochecker: bool = True
