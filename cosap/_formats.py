from dataclasses import dataclass


@dataclass
class FileFormats:
    TRIMMING_OUTPUT: str = "trimmed_{d[identification]}{d[read_no]}.fastq.gz"
    MAPPING_OUTPUT: str = "unprocessed_{identification}.bam"
    SORTING_OUTPUT: str = "sorted_{identification}.bam"
    INDEXING_OUTPUT: str = "indexed_{identification}.bam"
    MERGING_OUTPUT: str = "merged_{identification}.bam"
    MDUP_OUTPUT: str = "mdup_{identification}.bam"
    CALIBRATION_OUTPUT: str = "calibrated_{identification}.bam"
    CALIBRATION_TABLE: str = "calibration_table_{identification}.table"
    CALIBRATED_INDEXING_OUTPUT: str = "calibrated_indexed_{identification}.bam"
    GATK_UNFILTERED_OUTPUT: str = "all_{identification}.vcf"
    GATK_SNP_OUTPUT: str = "snp_{identification}.vcf"
    GATK_INDEL_OUTPUT: str = "indel_{identification}.vcf"
    GATK_OTHER_VARIANTS_OUTPUT: str = "other_variants_{identification}.vcf"
    ANNOTATING_OUTPUT: str = "annotated_{germline_identification}_{tumor_identification}.vcf"
