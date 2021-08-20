from dataclasses import dataclass


@dataclass
class FileFormats:
    MAPPING_OUTPUT: str = "{identification}_{algorithm}.bam"
    SORTING_OUTPUT: str = "sorted_{identification}_{algorithm}.bam"
    INDEXING_OUTPUT: str = "indexed_{identification}_{algorithm}.bam"
    MERGING_OUTPUT: str = "merged_{identification}_{algorithm}.bam"
    CALIBRATION_OUTPUT: str = "calibrated_{identification}_{algorithm}.bam"
    CALIBRATION_TABLE: str = "calibration_table_{identification}_{algorithm}.table"
    CALIBRATED_INDEXING_OUTPUT: str = (
        "calibrated_indexed_{identification}_{algorithm}.bam"
    )
    SAMTOOLS_PILEUP_OUTPUT: str = "pileup_{identification}_{algorithm}.vcf"
    GATK_SNP_OUTPUT: str = "snp_{identification}_{algorithm}.vcf"
    
