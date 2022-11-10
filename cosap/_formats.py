from dataclasses import dataclass

from cosap._utils import join_paths


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
    GATK_GVCF_OUTPUT: str = r"{identification}.g.vcf"
    ANNOTATING_OUTPUT: str = r"annotated_{identification}.vcf"
    ANNOVAR_OUTPUT: str = r"annotated_{identification}.avinput"
    QUALIMAP_PDF_OUTPUT: str = r"qualimap_{identification}.pdf"
    MOSDEPTH_OUTPUT: str = r"{identification}.mosdepth.summary.txt"


@dataclass
class FolderFormats:
    QUALIMAP_OUTPUT: str = r"qualimap_{identification}"

@dataclass
class OutputFolders:
    TRIMMING: str = "TRIMMED_FASTQ"
    MAPPING: str = "BAM"
    BAMQC: str = "BAMQC"
    PREPROCESSOR: str = "PREPROCESSOR"
    CALIBRATION: str = "CALIBRATED_BAM"
    VARIANT_CALLING: str = "VCF"
    ANNOTATION: str = "ANNOTATION"
    REPORT: str = "REPORT"
    LOG: str = "LOG"
    TEMP_OUTPUT: str = "TEMP"


@dataclass
class FolderedOutputs:
    GATK_SNP_OUTPUT: str = join_paths(
        OutputFolders.VARIANT_CALLING, "{library}", FileFormats.GATK_SNP_OUTPUT
    )
    MAPPING_OUTPUT: str = join_paths(
        OutputFolders.MAPPING, "{library}", FileFormats.MAPPING_OUTPUT
    )
    TRIMMING_OUTPUT: str = join_paths(
        OutputFolders.TRIMMING, FileFormats.TRIMMING_OUTPUT
    )
    MDUP_OUTPUT: str = join_paths(
        OutputFolders.PREPROCESSOR, "{library}", FileFormats.MDUP_OUTPUT
    )
    CALIBRATION_OUTPUT: str = join_paths(
        OutputFolders.CALIBRATION, FileFormats.CALIBRATION_OUTPUT
    )
    ELPREP_CALIBRATION_OUTPUT: str = join_paths(
        OutputFolders.CALIBRATION, FileFormats.ELPREP_CALIBRATION_OUTPUT
    )
    ANNOTATING_OUTPUT: str = join_paths(
        OutputFolders.ANNOTATION, "{library}", FileFormats.ANNOTATING_OUTPUT
    )
    
