from typing import Dict

import click

from .pipeline_runner import PipelineRunner
from .workflows import DNAPipeline, DNAPipelineInput


class Cosap:
    def run_pipeline(self, pipeline_config: Dict, snakemake: str):
        PipelineRunner.run_pipeline(pipeline_config=pipeline_config, backend=snakemake)

    def mock_run_pipeline(self, pipeline_config: Dict):
        # TODO: should this be here?
        raise NotImplementedError()


@click.group()
def cli():
    pass


@cli.command()
@cli.option(
    "--analysis_type",
    help="Type of the analysis, can be somatic or germline",
    required=True,
)
@cli.option("--workdir", help="Directory that outputs will be saved", required=True)
@cli.option("--normal_sample", help="Path to normal sample", required=False)
@cli.option(
    "--tumor_sample",
    help="Path to tumor sample. This option can be used multiple times.",
    required=False,
    multiple=True,
)
@cli.option("--bed_file", help="Path to BED file", required=False)
@cli.option(
    "--mapper",
    help="Mapper algorithm to use while aligning reads. This option can be used multiple times. Options = ['bwa', 'bwa2', 'bowtie2']",
    required=False,
    multiple=True,
)
@cli.option(
    "--variant_caller",
    help="Variant caller algorithm to use for variant detection. \
    This option can be used multiple times, Options = ['mutect2','varscan2','haplotypecaller','octopus','strelka','somaticsniper','vardict', 'deepvariant']",
    required=False,
    multiple=True,
)
@cli.option(
    "--normal_sample_name",
    help="Sample name of germline to write in bam. Default is 'normal'",
    required=False,
)
@cli.option(
    "--tumor_sample_name",
    help="Sample name of tumor to write in bam. Default is 'tumor'",
    required=False,
)
@cli.option(
    "--bam_qc",
    help="Qaulity control algorithm for .bam quality check. Options = ['qualimap', 'mosdepth']",
    required=False,
)
@cli.option(
    "--annotators",
    help="Annotation tool to annotate variants in vcf files.",
    required=False,
    multiple=True,
)
@cli.option(
    "--gvcf",
    help="Generate gvcf files",
    type=bool,
    required=False,
)
@cli.option(
    "--msi",
    help="Run microsatellite instability analysis",
    required=False,
    is_flag=True,
)
@cli.option(
    "--gene_fusion",
    help="Run gene fusion analysis",
    required=False,
    is_flag=True,
)
@cli.option(
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
