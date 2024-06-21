import glob
import os
from pathlib import Path

import yaml
from celery import Celery, shared_task

from .._config import AppConfig
from ..parsers import ProjectResultsParser
from ..pipeline_builder.builders import Annotator
from ..tools.annotators import AnnotatorFactory
from ..workflows import DNAPipeline, DNAPipelineInput
from ..workflows._variant_annotation import VariantMultipleAnnotator

celery_app = Celery("cosap")
celery_app.config_from_object("cosap.celery.celeryconfig")


@shared_task(name="dna_pipeline_task")
def dna_pipeline_task(
    analysis_type,
    workdir,
    normal_sample,
    tumor_samples,
    bed_file,
    mappers,
    variant_callers,
    annotation,
    bam_qc="qualimap",
    normal_sample_name="normal",
    tumor_sample_name="tumor",
    msi=True,
    gene_fusion=True,
):
    dna_pipeline_input = DNAPipelineInput(
        ANALYSIS_TYPE=analysis_type,
        WORKDIR=workdir,
        NORMAL_SAMPLE=normal_sample,
        TUMOR_SAMPLES=tumor_samples,
        BED_FILE=bed_file,
        MAPPERS=mappers,
        VARIANT_CALLERS=variant_callers,
        NORMAL_SAMPLE_NAME=normal_sample_name,
        TUMOR_SAMPLE_NAME=tumor_sample_name,
        BAM_QC=bam_qc,
        ANNOTATORS=annotation,
        MSI=msi,
        GENEFUSION=gene_fusion,
        DEVICE=AppConfig().DEVICE,
    )

    dna_pipeline = DNAPipeline(
        dna_pipeline_input=dna_pipeline_input,
    )
    return dna_pipeline.run_pipeline()


@shared_task(name="parse_project_results")
def parse_project_data(path):
    configs = glob.glob(f"{path}/*_config.yaml")
    latest_config = max(configs, key=os.path.getctime, default=None)
    if latest_config is None:
        raise FileNotFoundError
    config_dict = yaml.load(Path(latest_config).read_text(), Loader=yaml.Loader)
    parser = ProjectResultsParser(pipeline_config=config_dict)
    return {
        # "qc_coverage_histogram": parser.qc_coverage_histogram,
        "variants": parser.variants,
        "qc_results": parser.qc_genome_results,
        "msi_score": parser.msi_score,
    }


@shared_task(name="annotation_task")
def annotate_variants(variants: list, workdir: str) -> list:
    variant_annotator = VariantMultipleAnnotator(variants, workdir)
    return variant_annotator.annotate()
