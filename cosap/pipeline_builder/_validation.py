"""Pipeline configuration validation.

Validates pipeline config dicts before they are passed to runners,
catching misconfigurations early instead of deep inside tool execution.
"""

import os
import re
from typing import Dict, List, Optional


class PipelineConfigError(Exception):
    """Raised when a pipeline configuration is invalid."""

    pass


# Regex for safe identifiers: alphanumeric, underscore, hyphen, dot.
_SAFE_IDENTIFIER_RE = re.compile(r"^[a-zA-Z0-9._-]+$")

# Characters that must never appear in file paths used by the pipeline.
_DANGEROUS_PATH_CHARS = re.compile(r"[;\|`$\(\)\{\}!&<>]")

# Known valid library names per step type (lowercase).
VALID_MAPPERS = {"bwa", "bwa2", "bowtie", "bowtie2"}
VALID_VARIANT_CALLERS = {
    "mutect",
    "mutect2",
    "somaticsniper",
    "strelka",
    "varscan",
    "varscan2",
    "octopus",
    "muse",
    "vardict",
    "vardictjava",
    "haplotypecaller",
    "varscan_germline",
    "varnet",
    "manta",
    "deepvariant",
}
VALID_ANNOTATORS = {
    "annovar",
    "intervar",
    "cancervar",
    "vep",
    "pharmcat",
    "annotsv",
    "exomiser",
    "metarnn",
}
VALID_QC_TOOLS = {"qualimap", "mosdepth"}
VALID_PREPROCESSORS = {
    "mark_duplicate",
    "indexer",
    "base_recalibrator",
    "liftover_variant",
    "merger",
    "sorter",
    "trimmer",
    "elprep",
}
VALID_CNV_CALLERS = {"cnvkit", "gatk_germline_cnv"}
VALID_MSI_CALLERS = {"msisensor"}
VALID_GENE_FUSION_CALLERS = {"genefuse"}
VALID_MITOCHONDRIAL_CALLERS = {"gatk_mitochondrial"}
VALID_MODEL_TYPES = {"WGS", "WES"}


def validate_safe_path(path: str, field_name: str) -> None:
    """Validate that a path does not contain shell metacharacters or traversal."""
    if not path:
        return

    if _DANGEROUS_PATH_CHARS.search(path):
        raise PipelineConfigError(
            f"Invalid characters in {field_name}: '{path}'. "
            f"Paths must not contain shell metacharacters."
        )

    # Resolve to real path and check for traversal attempts.
    real_path = os.path.realpath(path)
    if ".." in path.split(os.sep):
        raise PipelineConfigError(
            f"Path traversal detected in {field_name}: '{path}'. "
            f"Use absolute paths without '..' components."
        )


def validate_identifier(value: str, field_name: str) -> None:
    """Validate that a string identifier is safe (no shell/path metacharacters)."""
    if not value:
        return

    if not _SAFE_IDENTIFIER_RE.match(value):
        raise PipelineConfigError(
            f"Invalid {field_name}: '{value}'. "
            f"Identifiers must contain only alphanumeric characters, "
            f"underscores, hyphens, and dots."
        )


def validate_library_name(library: str, valid_names: set, step_name: str) -> None:
    """Validate that a library name is recognized."""
    if not library:
        raise PipelineConfigError(f"Missing library name for {step_name} step.")

    if library.lower() not in valid_names:
        raise PipelineConfigError(
            f"Unknown {step_name} library: '{library}'. "
            f"Valid options are: {sorted(valid_names)}"
        )


def validate_file_exists(path: str, field_name: str) -> None:
    """Validate that a file path exists."""
    if not path:
        return

    resolved = os.path.realpath(path)
    if not os.path.exists(resolved):
        raise PipelineConfigError(
            f"File not found for {field_name}: '{path}' (resolved to '{resolved}')"
        )


def _validate_step_configs(
    step_configs: Dict,
    step_name: str,
    valid_libraries: Optional[set] = None,
) -> List[str]:
    """Validate all configs within a pipeline step. Returns list of errors."""
    errors = []

    if not isinstance(step_configs, dict):
        errors.append(
            f"{step_name} config must be a dictionary, got {type(step_configs).__name__}"
        )
        return errors

    for name, config in step_configs.items():
        # Validate the step name identifier
        try:
            validate_identifier(name, f"{step_name} step name")
        except PipelineConfigError as e:
            errors.append(str(e))

        if not isinstance(config, dict):
            continue

        # Validate library if applicable
        library = config.get("library")
        if library and valid_libraries:
            try:
                validate_library_name(library, valid_libraries, step_name)
            except PipelineConfigError as e:
                errors.append(str(e))

        # Validate file paths in config values
        path_keys = [
            "input",
            "output",
            "normal_input",
            "tumor_input",
            "input_bam",
            "input_bam_index",
            "bed",
            "table",
            "metrics",
            "reject_output",
        ]
        for key in path_keys:
            value = config.get(key)
            if isinstance(value, str) and value:
                try:
                    validate_safe_path(value, f"{step_name}.{name}.{key}")
                except PipelineConfigError as e:
                    errors.append(str(e))

        # Validate sample names
        sample_name_keys = ["germline_sample_name", "tumor_sample_name"]
        for key in sample_name_keys:
            value = config.get(key)
            if isinstance(value, str) and value:
                try:
                    validate_identifier(value, f"{step_name}.{name}.{key}")
                except PipelineConfigError as e:
                    errors.append(str(e))

        # Validate model_type for deepvariant
        library = config.get("library", "")
        if library.lower() == "deepvariant":
            model_type = config.get("model_type")
            if not model_type:
                errors.append(
                    f"{step_name}.{name}: model_type is required for deepvariant. "
                    f"Valid options are: {sorted(VALID_MODEL_TYPES)}"
                )
            elif model_type not in VALID_MODEL_TYPES:
                errors.append(
                    f"{step_name}.{name}: invalid model_type '{model_type}'. "
                    f"Valid options are: {sorted(VALID_MODEL_TYPES)}"
                )

    return errors


def validate_pipeline_config(pipeline_config: Dict) -> None:
    """Validate a complete pipeline configuration dict.

    Raises PipelineConfigError if any validation check fails.
    """
    errors = []

    if not isinstance(pipeline_config, dict):
        raise PipelineConfigError(
            f"Pipeline config must be a dictionary, got {type(pipeline_config).__name__}"
        )

    # Validate required top-level keys
    required_keys = ["workdir"]
    for key in required_keys:
        if key not in pipeline_config:
            errors.append(f"Missing required config key: '{key}'")

    # Validate workdir
    workdir = pipeline_config.get("workdir")
    if workdir:
        try:
            validate_safe_path(workdir, "workdir")
        except PipelineConfigError as e:
            errors.append(str(e))

    # Validate library_path
    library_path = pipeline_config.get("library_path")
    if library_path:
        try:
            validate_safe_path(library_path, "library_path")
        except PipelineConfigError as e:
            errors.append(str(e))

    # Validate each pipeline step section with its valid libraries
    step_validations = {
        "mapping": VALID_MAPPERS,
        "variant-calling": VALID_VARIANT_CALLERS,
        "annotation": VALID_ANNOTATORS,
        "quality_control": VALID_QC_TOOLS,
        "cnv": VALID_CNV_CALLERS,
        "germline_cnv": VALID_CNV_CALLERS,
        "mitochondrial_variant_calling": VALID_MITOCHONDRIAL_CALLERS,
    }

    for step_key, valid_libs in step_validations.items():
        step_config = pipeline_config.get(step_key, {})
        if step_config:
            errors.extend(_validate_step_configs(step_config, step_key, valid_libs))

    # Validate steps without library constraints
    for step_key in [
        "trim",
        "sorting",
        "index",
        "merge",
        "mdup",
        "calibrate",
        "elprep",
        "genefusion",
        "msi",
    ]:
        step_config = pipeline_config.get(step_key, {})
        if step_config:
            errors.extend(_validate_step_configs(step_config, step_key))

    if errors:
        error_msg = "Pipeline configuration validation failed:\n"
        for i, error in enumerate(errors, 1):
            error_msg += f"  {i}. {error}\n"
        raise PipelineConfigError(error_msg)
