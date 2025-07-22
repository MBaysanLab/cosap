# Tools
# Annotators
from ..tools.annotators._annotsv_annotator import AnnotSVAnnotator
from ..tools.annotators._annovar_annotator import AnnovarAnnotator
from ..tools.annotators._cancervar_annotatator import CancervarAnnotator
from ..tools.annotators._ensembl_vep_annotator import VepAnnotator
from ..tools.annotators._intervar_annotator import IntervarAnnotator
from ..tools.annotators._pharmcat_annotator import PharmcatAnnotator

# Special Callers
from ..tools.cnv_callers._cnvkit_cnv_caller import CNVKit
from ..tools.gene_fusion_callers._genefuse_fusion_caller import GeneFuse
from ..tools.msi_callers._msisensorpro_msicaller import MSISensorPro

# Mappers
from ..tools.mappers._bowtie_mapper import Bowtie2Mapper
from ..tools.mappers._bwa2_mapper import BWA2Mapper
from ..tools.mappers._bwa_mapper import BWAMapper

# Preprocessors
from ..tools.preprocessors._base_recalibrator import BaseRecalibrator
from ..tools.preprocessors._elprep_preprocess import ElprepPreprocess
from ..tools.preprocessors._indexer import BamIndexer
from ..tools.preprocessors._mark_duplicate import MarkDuplicate
from ..tools.preprocessors._sorter import SamtoolsSorter
from ..tools.preprocessors._split_bam_by_chr import SplitbyCHR
from ..tools.preprocessors._trimmer import Trimmer

# Quality Controllers
from ..tools.quality_controllers._mosdepth import Mosdepth
from ..tools.quality_controllers._qualimap import Qualimap

# Variant Callers
from ..tools.variant_callers._deepvariant_variantcaller import DeepVariantVariantCaller
from ..tools.variant_callers._haplotypecaller_variantcaller import HaplotypeCallerVariantCaller
from ..tools.variant_callers._manta_variantcaller import MantaVariantCaller
from ..tools.variant_callers._muse_variantcaller import MuseVariantCaller
from ..tools.variant_callers._mutect2_variantcaller import Mutect2VariantCaller
from ..tools.variant_callers._octopus_variantcaller import OctopusVariantCaller
from ..tools.variant_callers._somaticsniper_variantcaller import SomaticSniperVariantCaller
from ..tools.variant_callers._strelka_variantcaller import StrelkaVariantCaller
from ..tools.variant_callers._vardict_variantcaller import VarDictVariantCaller
from ..tools.variant_callers._varnet_variantcaller import VarNetVariantCaller
from ..tools.variant_callers._varscan_germline_variantcaller import VarScanGermlineVariantCaller
from ..tools.variant_callers._varscan_variantcaller import VarScanVariantCaller

from .downloadable_file import DownloadableFile


downloadable_files = [
    DownloadableFile(filename="Homo_sapiens_assembly38.dict",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict",
                     size=581712,
                     md5="3884c62eb0e53fa92459ed9bff133ae6",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta",
                     size=3249912778,
                     md5="7ff134953dcca8c8997453bbb80b6b5e",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.64.alt",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
                     size=487553,
                     md5="b07e65aa4425bc365141756f5c98328c",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.64.amb",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb",
                     size=20199,
                     md5="e4dc4fdb7358198e0847106599520aa9",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.64.ann",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann",
                     size=455474,
                     md5="af611ed0bb9487fb1ba4aa1a7e7ad21c",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.64.bwt",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt",
                     size=3217347004,
                     md5="7f0c8dcfc86b7c2ce3e3a54118d68fbd",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.64.pac",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac",
                     size=804336731,
                     md5="178862a79b043a2f974ef10e3877ef86",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.64.sa",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa",
                     size=1608673512,
                     md5="91a5d5ed3986db8a74782e5f4519eb5f",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.fai",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
                     size=160928,
                     md5="f76371b113734a56cde236bc0372de0a",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.dbsnp138.vcf",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
                     size=10950827213,
                     md5="f7e1ef5c1830bfb33675b9c7cbaa4868",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.dbsnp138.vcf.idx",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
                     size=12480412,
                     md5="2ba2e25b675c3a00ca663543b98170d8",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                     size=1888262073,
                     md5="b2979b47800b59b41920bf5432c4b2a0",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
                     size=2128536,
                     md5="1610d5349edae9df29bd0ce62fdea8ce",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                     size=20685880,
                     md5="2e02696032dcfe95ff0324f4a13508e3",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
                     url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
                     size=1500013,
                     md5="4c807e2cbe0752c0c44ac82ff3b52025",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.1.bt2",
                     url="https://drive.usercontent.google.com/download?id=1L4AQf6TsmoKPRyEEomf_tN3RkSzBpOwq&confirm=t",
                     size=1019133809,
                     md5="b0a021fcf1eedc19d43c43e2022f131e",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.2.bt2",
                     url="https://drive.usercontent.google.com/download?id=1PtCZReN1mfMvSW7rRcl1ZXXSrmEkp5Ln&confirm=t",
                     size=760863372,
                     md5="2f6632d01b0c9ee089e87537a2d7ab85",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.3.bt2",
                     url="https://drive.usercontent.google.com/download?id=1o5Z1kQ9hXXhkDxVW3S0T1StvVLpg1hSQ&confirm=t",
                     size=40751,
                     md5="131ab8d48cb890b85f52c726a9dd69a8",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.4.bt2",
                     url="https://drive.usercontent.google.com/download?id=1DLE1VWzCIj9t9HA4wlSWBajQFRjQqXMq&confirm=t",
                     size=760863367,
                     md5="c80305578802ebe1e06ff1f3c568c6e6",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.fasta.bwt.2bit.64",
                     url="https://drive.usercontent.google.com/download?id=1YL527ilA82bukLWSJ8Wr4zceoiQbWrqw&confirm=t",
                     size=10456377594,
                     md5="cbd28f3c971c1dd85ebce77ca471ed65",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.rev.1.bt2",
                     url="https://drive.usercontent.google.com/download?id=1Nbr6piXWOKTe1vrqryKBglLTtdXPr8Yo&confirm=t",
                     size=1019133809,
                     md5="e183ea30fc1d5412c0a8d266d3a41445",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.rev.2.bt2",
                     url="https://drive.usercontent.google.com/download?id=1x-oBuZBjOfkljg61VbXXr6hHLSwqpdTG&confirm=t",
                     size=760863372,
                     md5="bbe6bf20fd791075916cc0f07600e01d",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.dbsnp138.elsites",
                     url="https://drive.usercontent.google.com/download?id=1qO7YU5uLP3twmJFnq-3YrccIyovb7Af9&confirm=t",
                     size=1372365563,
                     md5="edd48d972df03967df20f52fc79eb0a5",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="Homo_sapiens_assembly38.elfasta",
                     url="https://drive.usercontent.google.com/download?id=1pqxlP1vErGZNCUIPaHcOP2iS5VmKsj37&confirm=t",
                     size=3217493596,
                     md5="e4c2886f35936aaef28f329b501858a0",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="genefuse_cancer.hg38.csv",
                     url="https://raw.githubusercontent.com/OpenGene/GeneFuse/master/genes/cancer.hg38.csv",
                     size=120406,
                     md5="8bda272f629bda6a8e44e79c24e59b9f",
                     requiring_steps={}
                     ),
    DownloadableFile(filename="msisensor_hg38.list",
                     url="https://raw.githubusercontent.com/xjtu-omics/msisensor-pro/3f8fa770f1cb329a49549750ac171f824c8dcadd/data/GRCh38.baseline_TCGA-v1.3.tsv",
                     size=889492,
                     md5="40d4e2f78b1fad4aa4e46fe4277473bc",
                     requiring_steps={}
                     ),
]
