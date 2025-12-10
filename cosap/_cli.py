import os

import click
from dotenv import set_key

from ._config import get_cosap_dotenv
from ._formats import FileFormats
from ._utils import join_paths
from .file_downloader import FileDownloader
from .pipeline_builder.builders import (
    MDUP,
    BamReader,
    FastqReader,
    GermlineCNVCaller,
    Mapper,
    Merger,
    Recalibrator,
    Sorter,
    Trimmer,
    VariantCaller,
)
from .tools import MapperFactory, PreprocessorFactory, VariantCallerFactory
from .tools.cnv_callers import CNVCallerFactory
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
@click.option(
    "--device",
    help="Device to run the mapping on. Options = ['cpu', 'gpu']",
    required=False,
    default="cpu",
)
def mapper(fastqs, mapper, output_dir, sample_name, read_groups, output_name, device):
    # Check if reads are provided
    if not fastqs:
        print("Please provide the path to the reads")
        exit()

    # Check if output_dir is provided
    if not output_dir:
        print("Current directory is used as the default output directory")
        output_dir = os.getcwd()

    if not sample_name and not read_groups:
        print("Please provide either a sample name or read groups")
        exit()

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
                "SM": sample_name,
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
    cosap_mapper.map(mapper_config[mapper_builder.key][mapper_builder.name], device=device)

    print(
        f"Output file: {join_paths(output_dir, FileFormats.MAPPING_OUTPUT.format(identification=output_name))}"
    )


@cli.command()
@click.option(
    "--variant_caller",
    help="Variant caller algorithm to use for variant detection. \
    Options = ['mutect2','varscan2','haplotypecaller','octopus','strelka','somaticsniper','vardict', 'deepvariant']",
    required=True,
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
    if tumor_sample:
        tumor_bam = BamReader(tumor_sample, name="tumor")
    if normal_sample:
        normal_bam = BamReader(normal_sample, name="normal")

    variant_caller_builder = VariantCaller(
        library=variant_caller,
        tumor_sample=tumor_bam if tumor_sample else None,
        normal_sample=normal_bam if normal_sample else None,
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


@cli.command()
@click.option(
    "--bam",
    "bam_files",
    help="Path to a BAM file. Can be specified multiple times for cohort analysis. "
         "Minimum 3 samples recommended for cohort mode.",
    required=True,
    multiple=True,
    type=click.Path(exists=True),
)
@click.option(
    "--intervals",
    help="Path to intervals file (BED or interval_list). Required for exome/panel data.",
    required=False,
    type=click.Path(exists=True),
)
@click.option(
    "--contig-ploidy-priors",
    "contig_ploidy_priors",
    help="Path to contig ploidy priors TSV file. Required for cohort mode.",
    required=False,
    type=click.Path(exists=True),
)
@click.option(
    "--output-dir",
    "output_dir",
    help="Output directory for results. Defaults to current directory.",
    required=False,
    type=click.Path(),
)
@click.option(
    "--name",
    help="Name for this analysis run. Used in output file naming.",
    required=False,
    default="germline_cnv_cohort",
)
@click.option(
    "--data-type",
    "data_type",
    help="Type of sequencing data.",
    type=click.Choice(["wgs", "exome"], case_sensitive=False),
    default="wgs",
)
@click.option(
    "--bin-length",
    "bin_length",
    help="Bin length for WGS data (default 1000). Set to 0 for exome data.",
    type=int,
    default=None,
)
@click.option(
    "--padding",
    help="Padding for intervals in bp (default 0 for WGS, 250 for exome).",
    type=int,
    default=None,
)
@click.option(
    "--sensitive",
    help="Use increased sensitivity parameters for WGS data. May increase false positives.",
    is_flag=True,
    default=False,
)
@click.option(
    "--allosomal-contig",
    "allosomal_contigs",
    help="Allosomal contig name (e.g., chrX, chrY). Can be specified multiple times.",
    multiple=True,
    default=["chrX", "chrY"],
)
@click.option(
    "--mappability-track",
    "mappability_track",
    help="Path to mappability track BED file for interval filtering.",
    required=False,
    type=click.Path(exists=True),
)
@click.option(
    "--max-copy-number",
    "max_copy_number",
    help="Maximum copy number to report (default 5).",
    type=int,
    default=5,
)
@click.option(
    "--verbosity",
    help="Logging verbosity level.",
    type=click.Choice(["INFO", "DEBUG", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
)
def germline_cnv_cohort(
    bam_files,
    intervals,
    contig_ploidy_priors,
    output_dir,
    name,
    data_type,
    bin_length,
    padding,
    sensitive,
    allosomal_contigs,
    mappability_track,
    max_copy_number,
    verbosity,
):
    """
    Run GATK GermlineCNVCaller pipeline in cohort mode.

    This command runs the complete GATK GermlineCNVCaller workflow for detecting
    rare germline copy number variants from multiple samples.

    The pipeline includes:

    \b
    Phase 1 - Cohort Model Building:
      1. PreprocessIntervals - Prepare intervals
      2. AnnotateIntervals - Annotate with GC content
      3. CollectReadCounts - Collect per-sample read counts
      4. FilterIntervals - Filter problematic intervals
      5. DetermineGermlineContigPloidy - Determine contig ploidies

    \b
    Phase 2 - Sample CNV Calling:
      6. GermlineCNVCaller - Call CNVs
      7. PostprocessGermlineCNVCalls - Generate per-sample VCFs

    \b
    Examples:
      # WGS cohort analysis
      cosap germline-cnv-cohort --bam sample1.bam --bam sample2.bam --bam sample3.bam \\
          --contig-ploidy-priors priors.tsv --output-dir ./results

      # Exome/panel analysis
      cosap germline-cnv-cohort --bam s1.bam --bam s2.bam --bam s3.bam \\
          --data-type exome --intervals targets.bed --output-dir ./results

      # WGS with increased sensitivity
      cosap germline-cnv-cohort --bam s1.bam --bam s2.bam --bam s3.bam \\
          --sensitive --output-dir ./results
    """
    # Validate inputs
    if len(bam_files) < 2:
        click.echo("Warning: Cohort mode works best with at least 3 samples. "
                   "Consider using more samples for better model training.", err=True)

    if data_type == "exome" and not intervals:
        click.echo("Error: --intervals is required for exome/panel data.", err=True)
        raise SystemExit(1)

    # Set output directory
    if not output_dir:
        output_dir = os.getcwd()
        click.echo(f"Using current directory as output: {output_dir}")
    
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Set default bin_length and padding based on data type
    if bin_length is None:
        bin_length = 0 if data_type == "exome" else 1000
    if padding is None:
        padding = 250 if data_type == "exome" else 0

    # Build params dict
    params = {
        "max_copy_number": max_copy_number,
        "verbosity": verbosity.upper(),
        "allosomal_contigs": list(allosomal_contigs),
    }

    if mappability_track:
        params["mappability_track"] = mappability_track

    # Add sensitive parameters if requested (for WGS)
    if sensitive:
        params.update({
            "class_coherence_length": 1000.0,
            "cnv_coherence_length": 1000.0,
            "enable_bias_factors": False,
            "interval_psi_scale": 1.0e-6,
            "log_mean_bias_standard_deviation": 0.01,
            "sample_psi_scale": 1.0e-6,
        })
        click.echo("Using increased sensitivity parameters for WGS data.")

    # For exome, ensure bias factors are enabled
    if data_type == "exome":
        params["enable_bias_factors"] = True

    # Create BAM readers
    click.echo(f"Processing {len(bam_files)} BAM files:")
    bam_readers = []
    for bam_path in bam_files:
        bam_name = os.path.basename(bam_path).replace(".bam", "")
        click.echo(f"  - {bam_name}: {bam_path}")
        bam_readers.append(BamReader(bam_path, name=bam_name))

    # Create GermlineCNVCaller builder
    germline_cnv_builder = GermlineCNVCaller(
        input_steps=bam_readers,
        intervals=intervals,
        contig_ploidy_priors=contig_ploidy_priors,
        name=name,
        bin_length=bin_length,
        padding=padding,
        allosomal_contigs=list(allosomal_contigs),
        params=params,
        output_dir=output_dir,
    )

    # Get config and run
    config = germline_cnv_builder.get_config()
    
    click.echo(f"\nStarting GATK GermlineCNVCaller pipeline...")
    click.echo(f"  Data type: {data_type.upper()}")
    click.echo(f"  Bin length: {bin_length}")
    click.echo(f"  Padding: {padding}")
    click.echo(f"  Output directory: {output_dir}")
    
    if intervals:
        click.echo(f"  Intervals: {intervals}")
    if contig_ploidy_priors:
        click.echo(f"  Contig ploidy priors: {contig_ploidy_priors}")

    # Create and run the CNV caller
    cnv_caller = CNVCallerFactory.create("gatk_germline_cnv")
    cnv_caller.call(config[germline_cnv_builder.key][germline_cnv_builder.name])

    click.echo(f"\nPipeline completed successfully!")
    click.echo(f"Results are in: {output_dir}")

