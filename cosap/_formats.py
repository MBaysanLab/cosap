from dataclasses import dataclass

from cosap._utils import join_paths


@dataclass
class FileFormats:
    TRIMMING_OUTPUT: str = r"trimmed_{identification}_{pair}.fastq.gz"
    TRIMMING_REPORT_OUTPUT: str = r"fastq_{identification}.json"
    MAPPING_OUTPUT: str = r"unprocessed_{identification}.bam"
    SORTING_OUTPUT: str = r"sorted_{identification}.bam"
    INDEXING_OUTPUT: str = r"{bam_file}.bai"
    MERGING_OUTPUT: str = r"merged_{identification}.bam"
    MDUP_OUTPUT: str = r"mdup_{identification}.bam"
    CALIBRATION_OUTPUT: str = r"calibrated_{identification}.bam"
    CALIBRATION_TABLE: str = r"calibration_table_{identification}.table"
    CALIBRATED_INDEXING_OUTPUT: str = r"calibrated_indexed_{identification}.bam"
    ELPREP_CALIBRATION_OUTPUT: str = r"elprep_calibrated_{identification}.bam"
    ALL_VARIANTS_OUTPUT: str = r"all_{identification}.vcf.gz"
    GATK_UNFILTERED_OUTPUT: str = r"all_unfiltered_{identification}.vcf.gz"
    SNP_OUTPUT: str = r"snp_{identification}.vcf.gz"
    INDEL_OUTPUT: str = r"indel_{identification}.vcf.gz"
    OTHER_VARIANTS_OUTPUT: str = r"other_variants_{identification}.vcf.gz"
    GVCF_OUTPUT: str = r"{identification}.g.vcf.gz"
    ANNOTATION_OUTPUT: str = r"annotated_{identification}.{custom_ext}"
    ANNOVAR_OUTPUT: str = r"annovar_{identification}.{sample}.avinput"
    QUALIMAP_OUTPUT: str = r"qualimap_{identification}"
    MOSDEPTH_OUTPUT: str = r"{identification}.mosdepth.summary.txt"
    SPLITTED_BAM_FILENAME: str = r"{name}_{split_no}.bam"
    GENEFUSION_OUTPUT: str = r"gene_fusion_{identification}.json"
    MSI_OUTPUT: str = r"msi_{identification}.txt"
    CNV_OUTPUT: str = r"cnv_{identification}.txt"


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
    GENE_FUSION: str = "GENE_FUSION"
    MSI: str = "MSI"
    CNV: str = "CNV"


@dataclass
class FolderedOutputs:
    VARIANT_CALLING_VCF_OUTPUT: str = join_paths(
        OutputFolders.VARIANT_CALLING, "{library}", FileFormats.ALL_VARIANTS_OUTPUT
    )
    VARIANT_CALLING_GVCF_OUTPUT: str = join_paths(
        OutputFolders.VARIANT_CALLING, "{library}", FileFormats.GVCF_OUTPUT
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
        OutputFolders.CALIBRATION, "{library}", FileFormats.ELPREP_CALIBRATION_OUTPUT
    )
    ANNOTATING_OUTPUT: str = join_paths(
        OutputFolders.ANNOTATION, "{library}", FileFormats.ANNOTATION_OUTPUT
    )
    BAMQC_OUTPUT: str = join_paths(OutputFolders.BAMQC, "{library}", "{identification}")
    REGIONS_FILE_OUTPUT: str = join_paths(OutputFolders.TEMP_OUTPUT, "regions")
    GENEFUSION_OUTPUT: str = join_paths(
        OutputFolders.GENE_FUSION, "{library}", FileFormats.GENEFUSION_OUTPUT
    )
    MSI_OUTPUT: str = join_paths(OutputFolders.MSI, "{library}", FileFormats.MSI_OUTPUT)
    CNV_OUTPUT: str = join_paths(OutputFolders.CNV, "{library}", FileFormats.CNV_OUTPUT)
