from typing import Dict

import click

from .default_pipelines import DNAPipeline
from .pipeline_runner import PipelineRunner


class Cosap:
    def run_pipeline(self, pipeline_config: Dict, snakemake: str):
        PipelineRunner.run_pipeline(pipeline_config=pipeline_config, backend=snakemake)

    def mock_run_pipeline(self, pipeline_config: Dict):
        # TODO: should this be here?
        raise NotImplementedError()


@click.command()
@click.option(
    "--analysis_type",
    help="Type of the analysis, can be somatic or germline",
    required=True,
)
@click.option("--workdir", help="Directory that outputs will be saved", required=True)
@click.option("--normal_sample", help="Path to normal sample", required=False)
@click.option(
    "--tumor_samples",
    help="Path to tumor sample. Can be multiple paths seperated with white space",
    required=False,
)
@click.option("--bed_file", help="Path to BED file", required=False)
@click.option(
    "--mappers",
    help="Mapper algorithms to use while aligning reads. Can be multiple option seperated with white space. Options = ['bwa', 'bwa2', 'bowtie2']",
    required=False,
)
@click.option(
    "--variant_callers",
    help="Variant caller algorithm to use for variant detection. \
    Can be multiple option seperated with white space, Options = ['mutect2','varscan2','haplotypecaller','octopus','strelka','somaticsniper','vardict', 'deepvariant']",
    required=False,
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
    help="Qaulity control algorithm for .bam quality check.",
    required=False,
)
@click.option(
    "--annotation",
    help="Annotation source to annotate variants in vcf files.",
    required=True,
)
def cosap_cli(
    analysis_type,
    workdir,
    normal_sample,
    tumor_samples,
    bed_file,
    mappers,
    variant_callers,
    normal_sample_name,
    tumor_sample_name,
    bam_qc,
    annotation,
):
    dna_pipeline = DNAPipeline(
        analysis_type=analysis_type,
        workdir=workdir,
        normal_sample=normal_sample,
        tumor_samples=tumor_samples,
        bed_file=bed_file,
        mappers=mappers,
        variant_callers=variant_callers,
        normal_sample_name=normal_sample_name,
        tumor_sample_name=tumor_sample_name,
        bam_qc=bam_qc,
        annotation=annotation,
    )
    dna_pipeline.run_pipeline()
