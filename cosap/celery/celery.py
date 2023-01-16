from celery import Celery, shared_task
from .._utils import parse_qualimap_coverage_histogram, parse_qualimap_genome_results
from ..default_pipelines import DNAPipeline

celery_app = Celery("cosap")
celery_app.config_from_object("cosap.celery.celeryconfig")

@shared_task(name="cosap_dna_pipeline_task")
def cosap_dna_pipeline_task(
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
    config_path = dna_pipeline.run_pipeline()
    return config_path

@shared_task(name="parse_qc_coverage_hist_task")
def parse_qc_coverage_hist_task(path):
    return parse_qualimap_coverage_histogram(path)


@shared_task(name="parse_qc_genome_results_task")
def parse_qc_genome_results_task(path):
    return parse_qualimap_genome_results(path)