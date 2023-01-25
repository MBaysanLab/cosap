from celery import Celery, shared_task
from ..parsers import ProjectResultsParser
from ..default_pipelines import DNAPipeline
import glob
import os
import yaml

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
    config = dna_pipeline.run_pipeline()
    return config

@shared_task(name="parse_project_results")
def parse_project_data(path):
    configs = glob.glob(f"{path}/*_config.yaml")
    latest_config = max(configs, key=os.path.getctime)

    config_dict = yaml.load(latest_config, Loader=yaml.Loader)
    parser = ProjectResultsParser(pipeline_config=config_dict)
    return {
        "qc_coverage_histogram": parser.qc_coverage_histogram,
        "variant_stats": parser.variant_stats,
        "variants": parser.variants,
        "qc_results": parser.qc_genome_results
    }