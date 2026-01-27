import glob
import json
import os
from logging import getLogger
from pathlib import Path
from tempfile import NamedTemporaryFile

import yaml
from celery import Celery, shared_task

from .._config import AppConfig
from .._utils import convert_vcf_to_json_fast
from ..parsers import ProjectResultsParser
from ..parsers.annotated_cnv_vcf_parser import AnnotatedCnvVcfParser
from ..workflows import (
    DNAPipeline,
    DNAPipelineInput,
    ExomiserPrioritizer,
    VariantMultipleAnnotator,
)

celery_app = Celery("cosap")
celery_app.config_from_object("cosap.celery.celeryconfig")

logger = getLogger(__name__)


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
    try:
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

        dna_pipeline.run_pipeline()
        return {"return_code": 0, "message": "Pipeline completed successfully"}

    except (FileNotFoundError, PermissionError, ValueError, Exception) as e:
        error_type = type(e).__name__
        logger.error(f"{error_type} in DNA pipeline: {e}")
        raise Exception(f"{error_type.replace('Error', '')}: {str(e)}")


@shared_task(name="parse_project_results")
def parse_project_data(path):
    try:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Project path does not exist: {path}")

        configs = glob.glob(f"{path}/*_config.yaml")
        latest_config = max(configs, key=os.path.getctime, default=None)

        if latest_config is None:
            raise FileNotFoundError(f"No config files found in {path}")

        config_dict = yaml.load(Path(latest_config).read_text(), Loader=yaml.Loader)

        parser = ProjectResultsParser(pipeline_config=config_dict)

        # Write variants to file
        with NamedTemporaryFile(mode="w", dir=path, suffix=".json", delete=False) as f:
            json.dump(json.loads(parser.variants), f)
            variants_file_path = f.name

        return {
            "variants": variants_file_path,
            "qc_results": parser.qc_genome_results,
            "msi_score": parser.msi_score,
        }

    except FileNotFoundError as e:
        logger.error(f"File not found error in parse_project_data: {e}")
        raise Exception(f"File not found: {str(e)}")
    except yaml.YAMLError as e:
        logger.error(f"YAML parsing error in parse_project_data: {e}")
        raise Exception(f"Invalid YAML configuration: {str(e)}")
    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error in parse_project_data: {e}")
        raise Exception(f"Invalid JSON data: {str(e)}")
    except PermissionError as e:
        logger.error(f"Permission error in parse_project_data: {e}")
        raise Exception(f"Permission denied: {str(e)}")
    except Exception as e:
        logger.error(f"Unexpected error in parse_project_data: {e}")
        raise Exception(f"Unexpected error: {str(e)}")


@shared_task(name="annotation_task")
def annotate_variants(
    variants_file_path: str, workdir: str, variant_type: str = None
) -> str:
    try:
        logger.info(
            f"Received request to annotate variants file: {variants_file_path}, workdir: {workdir}"
        )

        if not os.path.exists(variants_file_path):
            raise FileNotFoundError(f"Variants file not found: {variants_file_path}")

        if workdir is None:
            workdir = Path(variants_file_path).parent
        elif not os.path.exists(workdir):
            raise FileNotFoundError(f"Working directory not found: {workdir}")

        variant_annotator = VariantMultipleAnnotator(
            variants_file_path, workdir, variant_type=variant_type
        )
        variant_annotator.annotate()

        with NamedTemporaryFile(
            mode="w", dir=workdir, suffix=".json", delete=False
        ) as f:
            variant_annotator.save_results(f)
            annotated_variants_file_path = f.name

        return annotated_variants_file_path

    except FileNotFoundError as e:
        logger.error(f"File not found error in annotate_variants: {e}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error in annotate_variants: {e}")
        raise Exception(f"Invalid JSON in variants file: {str(e)}")
    except PermissionError as e:
        logger.error(f"Permission error in annotate_variants: {e}")
        raise Exception(f"Permission denied: {str(e)}")
    except Exception as e:
        logger.error(f"Unexpected error in annotate_variants: {e}")
        raise Exception(f"Annotation failed: {str(e)}")


@shared_task(name="parse_vcf_task")
def parse_vcf(
    vcf_path: str, sample_name: str, caller_type: str, reference_genome: str = None
) -> str:
    try:
        logger.info(
            f"Received request to parse VCF file: {vcf_path}, sample_name: {sample_name}, caller_type: {caller_type}"
        )

        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        # If reference is not hg38, we need to liftover the VCF
        if reference_genome and reference_genome.lower() != "hg38":
            from cosap.pipeline_builder import LiftoverVariant, VCFReader
            from cosap.tools.preprocessors import PreprocessorFactory

            liftover_builder = LiftoverVariant(
                input_step=VCFReader(vcf_path),
                name=f"{Path(vcf_path).stem}_liftover",
                output_dir=Path(vcf_path).parent,
            )
            liftover_preprocessor = PreprocessorFactory.create(
                preprocessor_type="liftover_variant",
            )
            liftover_preprocessor.run_preprocessor(
                liftover_config=liftover_builder.get_config(core_only=True),
            )
            vcf_path = liftover_builder.get_output()

        output_json = convert_vcf_to_json_fast(
            vcf_path, caller_type=caller_type, sample_name=sample_name
        )

        output_dir = Path(vcf_path).parent

        if not os.path.exists(output_dir):
            raise FileNotFoundError(f"Output directory not found: {output_dir}")

        # Write variants to file
        with NamedTemporaryFile(
            mode="w", dir=output_dir, suffix=".json", delete=False
        ) as f:
            json.dump(json.loads(output_json), f)
            variants_file_path = f.name

        return variants_file_path

    except FileNotFoundError as e:
        logger.error(f"File not found error in parse_vcf: {e}")
        raise
    except ValueError as e:
        logger.error(f"Invalid parameter in parse_vcf: {e}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error in parse_vcf: {e}")
        raise Exception(f"Invalid JSON output from VCF conversion: {str(e)}")
    except PermissionError as e:
        logger.error(f"Permission error in parse_vcf: {e}")
        raise Exception(f"Permission denied: {str(e)}")
    except Exception as e:
        logger.error(f"Unexpected error in parse_vcf: {e}")
        raise Exception(f"VCF parsing failed: {str(e)}")


@shared_task(name="exomiser_task")
def exomiser_task(vcf_path: str, hpo_terms: list, workdir: str) -> str:
    if workdir is None:
        workdir = str(Path(vcf_path).parent)

    variants_json = ExomiserPrioritizer(vcf_path, hpo_terms, workdir).prioritize()

    with NamedTemporaryFile(mode="w", dir=workdir, suffix=".json", delete=False) as f:
        json.dump(json.loads(variants_json), f)
        exomiser_variants_file_path = f.name

    return exomiser_variants_file_path


@shared_task(name="parse_annotated_cnv_task")
def parse_annotated_cnv_task(vcf_path: str, workdir: str = None) -> str:
    """
    Parse a VEP-annotated CNV VCF into a JSON list using AnnotatedCnvVcfParser.
    Returns the path to the written JSON file.
    """
    try:
        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        if workdir is None:
            workdir = str(Path(vcf_path).parent)
        elif not os.path.exists(workdir):
            raise FileNotFoundError(f"Working directory not found: {workdir}")

        parser = AnnotatedCnvVcfParser()
        parsed = parser.parse(vcf_path)

        with NamedTemporaryFile(mode="w", dir=workdir, suffix=".json", delete=False) as f:
            json.dump(parsed, f)
            return f.name

    except Exception as e:
        logger.error(f"parse_annotated_cnv_task failed: {e}")
        raise
