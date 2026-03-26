import glob
import json
import os
import re
from logging import getLogger
from pathlib import Path
from tempfile import NamedTemporaryFile

import yaml
from celery import Celery, shared_task
from celery.exceptions import SoftTimeLimitExceeded

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

# Regex for safe identifiers (sample names, caller types, etc.)
_SAFE_IDENTIFIER_RE = re.compile(r"^[a-zA-Z0-9._-]+$")

# Shell metacharacters that must never appear in paths
_DANGEROUS_PATH_CHARS = re.compile(r"[;\|`$\(\)\{\}!&<>]")


def _validate_path(path: str, field_name: str) -> str:
    """Validate and resolve a filesystem path.

    - Rejects paths with shell metacharacters.
    - Resolves to a real path (following symlinks).
    - Rejects paths with '..' components before resolution.

    Returns the resolved real path.
    Raises ValueError if the path is invalid.
    """
    if not path or not isinstance(path, str):
        raise ValueError(f"{field_name} must be a non-empty string, got: {path!r}")

    if _DANGEROUS_PATH_CHARS.search(path):
        raise ValueError(
            f"Invalid characters in {field_name}: path contains shell metacharacters."
        )

    # Reject explicit traversal attempts
    if ".." in path.split(os.sep):
        raise ValueError(
            f"Path traversal detected in {field_name}. "
            f"Use absolute paths without '..' components."
        )

    return os.path.realpath(path)


def _validate_identifier(value: str, field_name: str) -> str:
    """Validate that a string is a safe identifier."""
    if not value or not isinstance(value, str):
        raise ValueError(f"{field_name} must be a non-empty string, got: {value!r}")

    if not _SAFE_IDENTIFIER_RE.match(value):
        raise ValueError(
            f"Invalid {field_name}: '{value}'. "
            f"Must contain only alphanumeric characters, underscores, hyphens, and dots."
        )
    return value


def _validate_list_of_strings(values, field_name: str) -> list:
    """Validate that a value is a list of non-empty strings."""
    if not isinstance(values, list):
        raise ValueError(f"{field_name} must be a list, got: {type(values).__name__}")
    for i, v in enumerate(values):
        if not isinstance(v, str) or not v:
            raise ValueError(
                f"{field_name}[{i}] must be a non-empty string, got: {v!r}"
            )
    return values


@shared_task(
    bind=True,
    name="dna_pipeline_task",
    time_limit=172800,  # 48 hours hard limit
    soft_time_limit=169200,  # 47 hours soft limit
    acks_late=True,
    reject_on_worker_lost=True,
    max_retries=1,
    default_retry_delay=300,  # 5 minutes
)
def dna_pipeline_task(
    self,
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
        # Validate inputs at the task boundary
        _validate_identifier(analysis_type, "analysis_type")
        workdir = _validate_path(workdir, "workdir")
        _validate_identifier(normal_sample_name, "normal_sample_name")
        _validate_identifier(tumor_sample_name, "tumor_sample_name")

        if normal_sample:
            if isinstance(normal_sample, list):
                normal_sample = [
                    _validate_path(p, "normal_sample") for p in normal_sample
                ]
            elif isinstance(normal_sample, str):
                normal_sample = _validate_path(normal_sample, "normal_sample")

        if tumor_samples:
            if isinstance(tumor_samples, list):
                tumor_samples = [
                    _validate_path(p, "tumor_samples") for p in tumor_samples
                ]
            elif isinstance(tumor_samples, str):
                tumor_samples = _validate_path(tumor_samples, "tumor_samples")

        if bed_file:
            bed_file = _validate_path(bed_file, "bed_file")

        if mappers:
            _validate_list_of_strings(mappers, "mappers")
            for m in mappers:
                _validate_identifier(m, "mapper")

        if variant_callers:
            _validate_list_of_strings(variant_callers, "variant_callers")
            for vc in variant_callers:
                _validate_identifier(vc, "variant_caller")

        if bam_qc:
            _validate_identifier(bam_qc, "bam_qc")

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

    except SoftTimeLimitExceeded:
        logger.error("DNA pipeline task exceeded soft time limit")
        raise
    except (FileNotFoundError, PermissionError, ValueError, Exception) as e:
        error_type = type(e).__name__
        logger.error(f"{error_type} in DNA pipeline: {e}")
        raise Exception(f"{error_type.replace('Error', '')}: {str(e)}")


@shared_task(
    bind=True,
    name="parse_project_results",
    time_limit=3600,  # 1 hour hard limit
    soft_time_limit=3300,  # 55 minutes soft limit
    acks_late=True,
    reject_on_worker_lost=True,
    max_retries=2,
    default_retry_delay=60,
)
def parse_project_data(self, path):
    try:
        path = _validate_path(path, "path")

        if not os.path.exists(path):
            raise FileNotFoundError(f"Project path does not exist: {path}")

        if not os.path.isdir(path):
            raise ValueError(f"Project path is not a directory: {path}")

        configs = glob.glob(os.path.join(path, "*_config.yaml"))
        latest_config = max(configs, key=os.path.getctime, default=None)

        if latest_config is None:
            raise FileNotFoundError(f"No config files found in {path}")

        config_dict = yaml.safe_load(Path(latest_config).read_text())

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


@shared_task(
    bind=True,
    name="annotation_task",
    time_limit=43200,  # 12 hours hard limit
    soft_time_limit=41400,  # 11.5 hours soft limit
    acks_late=True,
    reject_on_worker_lost=True,
    max_retries=1,
    default_retry_delay=120,
)
def annotate_variants(
    self, variants_file_path: str, workdir: str, variant_type: str = None
) -> str:
    try:
        variants_file_path = _validate_path(variants_file_path, "variants_file_path")

        if variant_type is not None:
            _validate_identifier(variant_type, "variant_type")

        logger.info(
            f"Received request to annotate variants file: {variants_file_path}, workdir: {workdir}"
        )

        if not os.path.exists(variants_file_path):
            raise FileNotFoundError(f"Variants file not found: {variants_file_path}")

        if workdir is None:
            workdir = str(Path(variants_file_path).parent)
        else:
            workdir = _validate_path(workdir, "workdir")
            if not os.path.exists(workdir):
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


@shared_task(
    bind=True,
    name="parse_vcf_task",
    time_limit=3600,  # 1 hour hard limit
    soft_time_limit=3300,  # 55 minutes soft limit
    acks_late=True,
    reject_on_worker_lost=True,
    max_retries=2,
    default_retry_delay=60,
)
def parse_vcf(
    self,
    vcf_path: str,
    sample_name: str,
    caller_type: str,
    reference_genome: str = None,
) -> str:
    try:
        vcf_path = _validate_path(vcf_path, "vcf_path")
        _validate_identifier(sample_name, "sample_name")
        _validate_identifier(caller_type, "caller_type")
        if reference_genome is not None:
            _validate_identifier(reference_genome, "reference_genome")

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


@shared_task(
    bind=True,
    name="exomiser_task",
    time_limit=7200,  # 2 hours hard limit
    soft_time_limit=6900,  # 1 hour 55 minutes soft limit
    acks_late=True,
    reject_on_worker_lost=True,
    max_retries=1,
    default_retry_delay=120,
)
def exomiser_task(self, vcf_path: str, hpo_terms: list, workdir: str) -> str:
    vcf_path = _validate_path(vcf_path, "vcf_path")
    _validate_list_of_strings(hpo_terms, "hpo_terms")
    for term in hpo_terms:
        _validate_identifier(term, "hpo_term")

    if workdir is None:
        workdir = str(Path(vcf_path).parent)
    else:
        workdir = _validate_path(workdir, "workdir")

    variants_json = ExomiserPrioritizer(vcf_path, hpo_terms, workdir).prioritize()

    with NamedTemporaryFile(mode="w", dir=workdir, suffix=".json", delete=False) as f:
        json.dump(json.loads(variants_json), f)
        exomiser_variants_file_path = f.name

    return exomiser_variants_file_path


@shared_task(
    bind=True,
    name="parse_annotated_cnv_task",
    time_limit=3600,  # 1 hour hard limit
    soft_time_limit=3300,  # 55 minutes soft limit
    acks_late=True,
    reject_on_worker_lost=True,
    max_retries=2,
    default_retry_delay=60,
)
def parse_annotated_cnv_task(self, vcf_path: str, workdir: str = None) -> str:
    """
    Parse a VEP-annotated CNV VCF into a JSON list using AnnotatedCnvVcfParser.
    Returns the path to the written JSON file.
    """
    try:
        vcf_path = _validate_path(vcf_path, "vcf_path")

        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        if workdir is None:
            workdir = str(Path(vcf_path).parent)
        else:
            workdir = _validate_path(workdir, "workdir")
            if not os.path.exists(workdir):
                raise FileNotFoundError(f"Working directory not found: {workdir}")

        parser = AnnotatedCnvVcfParser()
        parsed = parser.parse(vcf_path)

        with NamedTemporaryFile(
            mode="w", dir=workdir, suffix=".json", delete=False
        ) as f:
            json.dump(parsed, f)
            return f.name

    except Exception as e:
        logger.error(f"parse_annotated_cnv_task failed: {e}")
        raise
