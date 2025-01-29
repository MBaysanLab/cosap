import os

import click
from dotenv import load_dotenv, set_key

from ._config import get_cosap_dotenv
from ._formats import FileFormats
from ._utils import join_paths
from .file_downloader import FileDownloader
from .pipeline_builder.builders import (MDUP, BamReader, FastqReader, Mapper,
                                        Merger, Recalibrator, Sorter, Trimmer,
                                        VariantCaller)
from .tools import MapperFactory, PreprocessorFactory, VariantCallerFactory
from .workflows import DNAPipeline, DNAPipelineInput


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--analysis_type",
    help="Type of the analysis, can be somatic or germline",
    required=True,
)
@click.option("--workdir", help="Directory that outputs will be saved", required=True)
@click.option("--normal_sample", help="Path to normal sample", required=False)
@click.option(
    "--tumor_sample",
    help="Path to tumor sample. This option can be used multiple times.",
    required=False,
    multiple=True,
)
@click.option("--bed_file", help="Path to BED file", required=False)
@click.option(
    "--mapper",
    help="Mapper algorithm to use while aligning reads. This option can be used multiple times. Options = ['bwa', 'bwa2', 'bowtie2']",
    required=False,
    multiple=True,
)
@click.option(
    "--variant_caller",
    help="Variant caller algorithm to use for variant detection. \
    This option can be used multiple times, Options = ['mutect2','varscan2','haplotypecaller','octopus','strelka','somaticsniper','vardict', 'deepvariant']",
    required=False,
    multiple=True,
)
@click.option(
    "--normal_sample_name",
    help="Sample name of germline to write in bam. Default is 'normal'",
    required=False,
)
@click.option(
    "--tumor_sample_name",
    help="Sample name of tumor to write in bam. Default is 'tumor'",
    required=False,
)
@click.option(
    "--bam_qc",
    help="Qaulity control algorithm for .bam quality check. Options = ['qualimap', 'mosdepth']",
    required=False,
)
@click.option(
    "--annotators",
    help="Annotation tool to annotate variants in vcf files.",
    required=False,
    multiple=True,
)
@click.option(
    "--gvcf",
    help="Generate gvcf files",
    type=bool,
    required=False,
)
@click.option(
    "--msi",
    help="Run microsatellite instability analysis",
    required=False,
    is_flag=True,
)
@click.option(
    "--gene_fusion",
    help="Run gene fusion analysis",
    required=False,
    is_flag=True,
)
@click.option(
    "--device",
    help="Device to run the pipeline on. Options = ['cpu', 'gpu']",
    required=False,
    default="cpu",
)
def pipeline(
    analysis_type,
    workdir,
    normal_sample,
    tumor_sample,
    bed_file,
    mapper,
    variant_caller,
    normal_sample_name,
    tumor_sample_name,
    bam_qc,
    annotators,
    gvcf,
    msi,
    gene_fusion,
    device=None,
):

    print(
        f"Running cosap pipeline with the following parameters: \n \
            analysis_type: {analysis_type} \n \
            workdir: {workdir} \n \
            normal_sample: {normal_sample} \n \
            tumor_samples: {tumor_sample} \n \
            bed_file: {bed_file} \n \
            mappers: {mapper} \n \
            variant_callers: {variant_caller} \n \
            normal_sample_name: {normal_sample_name} \n \
            tumor_sample_name: {tumor_sample_name} \n \
            bam_qc: {bam_qc} \n \
            annotators: {annotators} \n \
            gvcf: {gvcf} \n \
            msi: {msi} \n \
            gene_fusion: {gene_fusion} \n \
            device: {device} \n"
    )

    pipeline_input = DNAPipelineInput(
        ANALYSIS_TYPE=analysis_type,
        WORKDIR=workdir,
        NORMAL_SAMPLE=tuple(normal_sample.split(",")) if normal_sample else None,
        TUMOR_SAMPLES=[tuple(t_sm.split(",")) for t_sm in tumor_sample],
        BED_FILE=bed_file,
        MAPPERS=mapper,
        VARIANT_CALLERS=variant_caller,
        NORMAL_SAMPLE_NAME=normal_sample_name if normal_sample_name else None,
        TUMOR_SAMPLE_NAME=tumor_sample_name if tumor_sample_name else None,
        BAM_QC=bam_qc,
        ANNOTATORS=annotators,
        GVCF=gvcf,
        MSI=msi,
        GENEFUSION=gene_fusion,
        DEVICE=device,
    )
    dna_pipeline = DNAPipeline(pipeline_input)
    dna_pipeline.run_pipeline()


@cli.command()
def download_data():
    # Check if COSAP_LIBRARY_PATH is set
    if "COSAP_LIBRARY_PATH" not in os.environ:
        print(
            "COSAP_LIBRARY_PATH is not set. Please set the path to the reference genome directory using `cosap set-ref-dir`"
        )
        exit()

    # Download the data
    file_downloader = FileDownloader()
    file_downloader.download_files()


@cli.command()
@click.option(
    "--path",
    help="Path to the reference genome directory",
    required=True,
)
def set_ref_dir(path):

    if existing_path := os.environ.get("COSAP_LIBRARY_PATH"):
        if input(
            f"This environment variable is already set to {existing_path}. Continue? [Y/n]: "
        ) not in [
            "",
            "y",
            "Y",
        ]:
            print("Exiting...")
            exit()

    # Set the path
    set_key(get_cosap_dotenv(), "COSAP_LIBRARY_PATH", path)


@cli.command()
@click.option(
    "--fastqs",
    help="White space separated paths to the reads",
    required=True,
    nargs=2,
)
@click.option(
    "--mapper",
    help="Mapper algorithm to use while aligning reads. Options = ['bwa', 'bwa2', 'bowtie2']",
    required=True,
)
@click.option(
    "--output_dir",
    help="Output directory to save the mapped reads",
    required=False,
)
@click.option(
    "--sample_name",
    help="Name of the sample",
    required=False,
)
@click.option(
    "--read_groups",
    help="Read groups to add to the bam file",
    required=False,
    multiple=True,
)
@click.option(
    "--output_name",
    help="Name of the output file",
    required=False,
)
def mapper(fastqs, mapper, output_dir, sample_name, read_groups, output_name):

    # Check if reads are provided
    if not fastqs:
        print("Please provide the path to the reads")
        exit()

    # Check if output_dir is provided
    if not output_dir:
        print("Current directory is used as the default output directory")
        output_dir = os.getcwd()

    fastqs = [
        FastqReader(fastq, name=sample_name, read=i + 1)
        for i, fastq in enumerate(fastqs)
    ]

    # Parse read groups as a dictionary
    read_groups = dict([rg.split(":") for rg in read_groups])
    # Run the pipeline
    mapper_builder = Mapper(
        input_step=fastqs,
        library=mapper,
        params={
            "read_groups": {
                "ID": "0",
                "SM": "test",
                "PU": "0",
                "PL": "ILLUMINA",
                "LB": "0",
            }
        },
        output_dir=output_dir,
        name=output_name,
    )
    mapper_config = mapper_builder.get_config()
    cosap_mapper = MapperFactory.create(mapper)
    cosap_mapper.map(mapper_config[mapper_builder.key][mapper_builder.name])

    print(
        f"Output file: {join_paths(output_dir,FileFormats.MAPPING_OUTPUT.format(identification=output_name))}"
    )


@cli.command()
@click.option(
    "--variant_caller",
    help="Variant caller algorithm to use for variant detection. \
    Options = ['mutect2','varscan2','haplotypecaller','octopus','strelka','somaticsniper','vardict', 'deepvariant']",
    required=True,
    multiple=True,
)
@click.option(
    "--tumor_sample",
    help="Path to tumor sample",
    required=False,
)
@click.option(
    "--normal_sample",
    help="Path to normal sample",
    required=False,
)
@click.option(
    "--output_dir",
    help="Output directory to save the variants",
    required=False,
)
@click.option(
    "--tumor_sample_name",
    help="Name of the tumor sample",
    required=False,
)
def variant_caller(
    variant_caller,
    tumor_sample,
    normal_sample,
    output_dir,
    tumor_sample_name,
):

    # Check if tumor_sample is provided
    if not tumor_sample and not normal_sample:
        print("Please provide the path to the analysis ready bam files.")
        exit()

    # Check if output_dir is provided
    if not output_dir:
        print("Current directory is used as the default output directory")
        output_dir = os.getcwd()

    # Create variant caller object with builder pattern
    tumor_bam = BamReader(tumor_sample, name="tumor")
    normal_bam = BamReader(normal_sample, name="normal")

    variant_caller_builder = VariantCaller(
        library=variant_caller,
        tumor_sample=tumor_bam,
        normal_sample=normal_bam,
        params={
            "tumor_sample_name": tumor_sample_name,
        },
        output_dir=output_dir,
    )
    variant_caller_config = variant_caller_builder.get_config()
    cosap_variant_caller = VariantCallerFactory.create(variant_caller)
    cosap_variant_caller.call_variants(
        variant_caller_config[variant_caller_builder.key][variant_caller_builder.name]
    )


@cli.command()
@click.option(
    "--input_file",
    help="Path to the input file",
    required=True,
)
@click.option(
    "--output_dir",
    help="Output directory to save the annotated file",
    required=False,
)
@click.option(
    "--preprocessor",
    help="Preprocessor bam files for variant calling. Options = ['mark_duplicate', 'base_recalibrator', 'merger', 'sorter', 'trimmer']",
    required=True,
)
def bam_preprocessor(
    preprocessor,
    input_file,
    output_dir,
):

    # Check if input_file is provided
    if not input_file:
        print("Please provide the path to the input file")
        exit()

    # Check if output_dir is provided
    if not output_dir:
        print("Current directory is used as the default output directory")
        output_dir = os.getcwd()

    input_file = BamReader(input_file)

    # Create preprocessor object with builder pattern
    def get_preprocessor(preprocessor):
        if preprocessor == "mark_duplicate":
            return MDUP(input_file, output_dir)
        elif preprocessor == "base_recalibrator":
            return Recalibrator(input_file, output_dir)
        elif preprocessor == "merger":
            return Merger(input_file, output_dir)
        elif preprocessor == "sorter":
            return Sorter(input_file, output_dir)
        elif preprocessor == "trimmer":
            return Trimmer(input_file, output_dir)
        else:
            print("Invalid preprocessor")
            exit()

    preprocessor_obj = get_preprocessor(preprocessor)
    proprocessor_config = preprocessor_obj.get_config()

    cosap_preprocessor = PreprocessorFactory.create(preprocessor)
    cosap_preprocessor.run_preprocessor(
        proprocessor_config[preprocessor_obj.key][preprocessor_obj.name]
    )
