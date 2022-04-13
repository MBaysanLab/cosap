from dataclasses import dataclass


@dataclass
class FileFormats:
    TRIMMING_OUTPUT: str = r"trimmed_{identification}_{pair}.fastq.gz"
    MAPPING_OUTPUT: str = r"unprocessed_{identification}.bam"
    SORTING_OUTPUT: str = r"sorted_{identification}.bam"
    INDEXING_OUTPUT: str = r"{prefix}_{identification}.bai"
    MERGING_OUTPUT: str = r"merged_{identification}.bam"
    MDUP_OUTPUT: str = r"mdup_{identification}.bam"
    CALIBRATION_OUTPUT: str = r"calibrated_{identification}.bam"
    CALIBRATION_TABLE: str = r"calibration_table_{identification}.table"
    CALIBRATED_INDEXING_OUTPUT: str = r"calibrated_indexed_{identification}.bam"
    ELPREP_CALIBRATION_OUTPUT: str = r"elprep_calibrated_{identification}.bam"
    GATK_UNFILTERED_OUTPUT: str = r"all_{identification}.vcf"
    GATK_FILTERED_OUTPUT: str = r"all_w_filters_{identification}.vcf"
    GATK_SNP_OUTPUT: str = r"snp_{identification}.vcf"
    GATK_INDEL_OUTPUT: str = r"indel_{identification}.vcf"
    GATK_OTHER_VARIANTS_OUTPUT: str = r"other_variants_{identification}.vcf"
    ANNOTATING_OUTPUT: str = r"annotated_{identification}.vcf"
    ANNOVAR_OUTPUT: str = r"annotated_{identification}.avinput"
