from dataclasses import dataclass


@dataclass
class FileFormats:
    MAPPING_OUTPUT: str = "{identification}.bam"
    SORTING_OUTPUT: str = "sorted_{identification}.bam"
    INDEXING_OUTPUT: str = "indexed_{identification}.bam"
    MERGING_OUTPUT: str = "merged_{identification}.bam"
    CALIBRATION_OUTPUT: str = "calibrated_{identification}.bam"
    CALIBRATION_TABLE: str = "calibration_table_{identification}.table"
    CALIBRATED_INDEXING_OUTPUT: str = "calibrated_indexed_{identification}.bam"
    SAMTOOLS_PILEUP_OUTPUT: str = "pileup_{identification}.vcf"
    GATK_SNP_OUTPUT: str = "snp_{identification}.vcf"
