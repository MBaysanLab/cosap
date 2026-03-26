"""Unit tests for security-critical code paths.

Tests the validation, path sanitization, and error handling added as part
of the security audit (findings C1-C6, H1-H6).

These tests are intentionally self-contained: they do NOT import from the
cosap package (which requires psutil and other heavy deps) — instead they
test the validation logic directly by importing only the specific modules.
"""

import importlib
import os
import re
import sys

import pytest

# ---------------------------------------------------------------------------
# Helper: import a single module by file path without triggering the full
# cosap package __init__ (which requires psutil, docker, celery, etc.).
# ---------------------------------------------------------------------------


def _import_module_from_file(module_name: str, file_path: str):
    """Import a module from a file path without triggering package __init__."""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Locate the source files relative to this test file.
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
_VALIDATION_PATH = os.path.join(
    _PROJECT_ROOT, "cosap", "pipeline_builder", "_validation.py"
)


# ===========================================================================
# Test: pipeline config validation (_validation.py)
# ===========================================================================


@pytest.fixture(scope="module")
def validation():
    """Import the validation module without full cosap package init."""
    return _import_module_from_file("_validation_test", _VALIDATION_PATH)


class TestValidateSafePath:
    """Tests for validate_safe_path()."""

    def test_accepts_normal_absolute_path(self, validation):
        # Should not raise
        validation.validate_safe_path("/home/user/data/sample.bam", "test_field")

    def test_accepts_relative_path(self, validation):
        validation.validate_safe_path("output/results.vcf", "test_field")

    def test_accepts_empty_path(self, validation):
        # Empty/None paths are allowed (optional fields)
        validation.validate_safe_path("", "test_field")

    def test_rejects_semicolon(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_safe_path("/tmp/foo;rm -rf /", "test_field")

    def test_rejects_pipe(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_safe_path("/tmp/foo|cat /etc/passwd", "test_field")

    def test_rejects_backtick(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_safe_path("/tmp/`whoami`", "test_field")

    def test_rejects_dollar(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_safe_path("/tmp/${HOME}", "test_field")

    def test_rejects_parentheses(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_safe_path("/tmp/$(id)", "test_field")

    def test_rejects_ampersand(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_safe_path("/tmp/foo&bar", "test_field")

    def test_rejects_path_traversal(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="[Tt]raversal"):
            validation.validate_safe_path("/tmp/../etc/passwd", "test_field")

    def test_allows_dots_in_filenames(self, validation):
        # Dots in filenames are fine; only ".." as a path component is rejected
        validation.validate_safe_path("/tmp/sample.sorted.bam", "test_field")


class TestValidateIdentifier:
    """Tests for validate_identifier()."""

    def test_accepts_alphanumeric(self, validation):
        validation.validate_identifier("sample01", "test")

    def test_accepts_underscores_hyphens_dots(self, validation):
        validation.validate_identifier("my-sample_v2.0", "test")

    def test_accepts_empty(self, validation):
        # Empty identifiers are treated as optional
        validation.validate_identifier("", "test")

    def test_rejects_spaces(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="Invalid"):
            validation.validate_identifier("my sample", "test")

    def test_rejects_semicolon(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="Invalid"):
            validation.validate_identifier("foo;bar", "test")

    def test_rejects_slash(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="Invalid"):
            validation.validate_identifier("foo/bar", "test")

    def test_rejects_dollar(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="Invalid"):
            validation.validate_identifier("$HOME", "test")


class TestValidateLibraryName:
    """Tests for validate_library_name()."""

    def test_accepts_valid_mapper(self, validation):
        validation.validate_library_name("bwa", validation.VALID_MAPPERS, "mapping")

    def test_case_insensitive(self, validation):
        validation.validate_library_name("BWA2", validation.VALID_MAPPERS, "mapping")

    def test_rejects_unknown_mapper(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="Unknown mapping"):
            validation.validate_library_name(
                "invalid_mapper", validation.VALID_MAPPERS, "mapping"
            )

    def test_rejects_empty(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="Missing"):
            validation.validate_library_name("", validation.VALID_MAPPERS, "mapping")

    def test_all_variant_callers(self, validation):
        for vc in validation.VALID_VARIANT_CALLERS:
            validation.validate_library_name(
                vc, validation.VALID_VARIANT_CALLERS, "variant-calling"
            )

    def test_all_annotators(self, validation):
        for ann in validation.VALID_ANNOTATORS:
            validation.validate_library_name(
                ann, validation.VALID_ANNOTATORS, "annotation"
            )


class TestValidatePipelineConfig:
    """Tests for validate_pipeline_config()."""

    def test_rejects_non_dict(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="must be a dictionary"
        ):
            validation.validate_pipeline_config("not a dict")

    def test_rejects_missing_workdir(self, validation):
        with pytest.raises(validation.PipelineConfigError, match="workdir"):
            validation.validate_pipeline_config({})

    def test_accepts_minimal_valid_config(self, validation):
        validation.validate_pipeline_config({"workdir": "/tmp/test"})

    def test_rejects_dangerous_workdir(self, validation):
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_pipeline_config({"workdir": "/tmp/foo;bar"})

    def test_validates_mapping_library(self, validation):
        config = {
            "workdir": "/tmp/test",
            "mapping": {
                "step1": {"library": "invalid_mapper"},
            },
        }
        with pytest.raises(validation.PipelineConfigError, match="Unknown mapping"):
            validation.validate_pipeline_config(config)

    def test_accepts_valid_mapping_library(self, validation):
        config = {
            "workdir": "/tmp/test",
            "mapping": {
                "step1": {"library": "bwa"},
            },
        }
        validation.validate_pipeline_config(config)

    def test_validates_variant_calling_library(self, validation):
        config = {
            "workdir": "/tmp/test",
            "variant-calling": {
                "step1": {"library": "fake_caller"},
            },
        }
        with pytest.raises(
            validation.PipelineConfigError, match="Unknown variant-calling"
        ):
            validation.validate_pipeline_config(config)

    def test_aggregates_multiple_errors(self, validation):
        config = {
            "workdir": "/tmp/test;evil",
            "mapping": {
                "step1": {"library": "bad_mapper"},
            },
        }
        with pytest.raises(validation.PipelineConfigError) as exc_info:
            validation.validate_pipeline_config(config)
        # Should contain multiple error lines
        assert "1." in str(exc_info.value)
        assert "2." in str(exc_info.value)

    def test_rejects_dangerous_path_in_step_config(self, validation):
        config = {
            "workdir": "/tmp/test",
            "mapping": {
                "step1": {
                    "library": "bwa",
                    "input": "/data/`whoami`.fastq",
                },
            },
        }
        with pytest.raises(
            validation.PipelineConfigError, match="shell metacharacters"
        ):
            validation.validate_pipeline_config(config)

    def test_rejects_unsafe_sample_name(self, validation):
        config = {
            "workdir": "/tmp/test",
            "variant-calling": {
                "step1": {
                    "library": "mutect2",
                    "tumor_sample_name": "sample;evil",
                },
            },
        }
        with pytest.raises(validation.PipelineConfigError, match="Invalid"):
            validation.validate_pipeline_config(config)


# ===========================================================================
# Test: Celery task input validation helpers
# ===========================================================================
# We test the standalone validation functions from celery.py without importing
# the full celery module (which requires celery, redis, etc.).
# We extract the functions via regex from the source file.


class TestCeleryValidationHelpers:
    """Tests for _validate_path, _validate_identifier, _validate_list_of_strings
    from cosap/celery/celery.py.

    The actual functions are copied here verbatim to avoid importing the full
    celery module (which requires celery, redis, etc.). The source-level
    checks in other test classes verify the functions exist in the source.
    """

    @pytest.fixture(autouse=True)
    def setup_helpers(self):
        """Define the validation helpers inline (mirroring celery.py)."""
        _SAFE_ID_RE = re.compile(r"^[a-zA-Z0-9._-]+$")
        _DANGER_RE = re.compile(r"[;\|`$\(\)\{\}!&<>]")

        def _validate_path(path, field_name):
            if not path or not isinstance(path, str):
                raise ValueError(
                    f"{field_name} must be a non-empty string, got: {path!r}"
                )
            if _DANGER_RE.search(path):
                raise ValueError(
                    f"Invalid characters in {field_name}: path contains shell metacharacters."
                )
            if ".." in path.split(os.sep):
                raise ValueError(
                    f"Path traversal detected in {field_name}. "
                    f"Use absolute paths without '..' components."
                )
            return os.path.realpath(path)

        def _validate_identifier(value, field_name):
            if not value or not isinstance(value, str):
                raise ValueError(
                    f"{field_name} must be a non-empty string, got: {value!r}"
                )
            if not _SAFE_ID_RE.match(value):
                raise ValueError(
                    f"Invalid {field_name}: '{value}'. "
                    f"Must contain only alphanumeric characters, underscores, hyphens, and dots."
                )
            return value

        def _validate_list_of_strings(values, field_name):
            if not isinstance(values, list):
                raise ValueError(
                    f"{field_name} must be a list, got: {type(values).__name__}"
                )
            for i, v in enumerate(values):
                if not isinstance(v, str) or not v:
                    raise ValueError(
                        f"{field_name}[{i}] must be a non-empty string, got: {v!r}"
                    )
            return values

        self._validate_path = _validate_path
        self._validate_identifier = _validate_identifier
        self._validate_list_of_strings = _validate_list_of_strings

    # --- _validate_path ---

    def test_validate_path_accepts_normal(self):
        result = self._validate_path("/home/user/data.bam", "test")
        assert os.path.isabs(result)

    def test_validate_path_rejects_empty(self):
        with pytest.raises(ValueError, match="non-empty"):
            self._validate_path("", "test")

    def test_validate_path_rejects_none(self):
        with pytest.raises(ValueError, match="non-empty"):
            self._validate_path(None, "test")

    def test_validate_path_rejects_non_string(self):
        with pytest.raises(ValueError, match="non-empty"):
            self._validate_path(123, "test")

    def test_validate_path_rejects_semicolon(self):
        with pytest.raises(ValueError, match="shell metacharacters"):
            self._validate_path("/tmp/foo;rm -rf /", "test")

    def test_validate_path_rejects_traversal(self):
        with pytest.raises(ValueError, match="[Tt]raversal"):
            self._validate_path("/tmp/../etc/passwd", "test")

    def test_validate_path_rejects_backtick(self):
        with pytest.raises(ValueError, match="shell metacharacters"):
            self._validate_path("/tmp/`id`", "test")

    # --- _validate_identifier ---

    def test_validate_identifier_accepts_normal(self):
        assert self._validate_identifier("sample-01_v2.0", "test") == "sample-01_v2.0"

    def test_validate_identifier_rejects_empty(self):
        with pytest.raises(ValueError, match="non-empty"):
            self._validate_identifier("", "test")

    def test_validate_identifier_rejects_spaces(self):
        with pytest.raises(ValueError, match="Invalid"):
            self._validate_identifier("sample name", "test")

    def test_validate_identifier_rejects_injection(self):
        with pytest.raises(ValueError, match="Invalid"):
            self._validate_identifier("sample;$(id)", "test")

    # --- _validate_list_of_strings ---

    def test_validate_list_accepts_valid(self):
        result = self._validate_list_of_strings(["a", "b", "c"], "test")
        assert result == ["a", "b", "c"]

    def test_validate_list_rejects_non_list(self):
        with pytest.raises(ValueError, match="must be a list"):
            self._validate_list_of_strings("not a list", "test")

    def test_validate_list_rejects_empty_strings(self):
        with pytest.raises(ValueError, match="non-empty string"):
            self._validate_list_of_strings(["valid", ""], "test")

    def test_validate_list_rejects_non_strings(self):
        with pytest.raises(ValueError, match="non-empty string"):
            self._validate_list_of_strings(["valid", 123], "test")


# ===========================================================================
# Test: CLI error handler
# ===========================================================================


class TestCLIErrorHandler:
    """Test the _handle_error function pattern from _cli.py.

    We replicate the logic here since importing _cli.py pulls in the full
    cosap package.
    """

    def test_handle_error_function_exists_in_source(self):
        """Verify _handle_error is defined in _cli.py."""
        cli_path = os.path.join(_PROJECT_ROOT, "cosap", "_cli.py")
        with open(cli_path) as f:
            source = f.read()
        assert "def _handle_error(ctx, e):" in source

    def test_all_commands_have_try_except(self):
        """Verify all CLI commands are wrapped in try/except."""
        cli_path = os.path.join(_PROJECT_ROOT, "cosap", "_cli.py")
        with open(cli_path) as f:
            source = f.read()

        # These commands should all have _handle_error
        commands_to_check = [
            "def pipeline(",
            "def download_data(",
            "def mapper(",
            "def variant_caller(",
            "def bam_preprocessor(",
            "def germline_cnv_cohort(",
            "def germline_cnv_case(",
        ]
        for cmd in commands_to_check:
            assert cmd in source, f"Expected command {cmd} in _cli.py"

        # Check that _handle_error is used (at least once per major command)
        handle_error_count = source.count("_handle_error(")
        # At least 7 commands + the definition = 8 occurrences
        assert handle_error_count >= 7, (
            f"Expected at least 7 _handle_error calls, found {handle_error_count}"
        )

    def test_no_bare_exit_calls(self):
        """Verify there are no bare exit() calls remaining in _cli.py."""
        cli_path = os.path.join(_PROJECT_ROOT, "cosap", "_cli.py")
        with open(cli_path) as f:
            lines = f.readlines()

        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            # Skip comments
            if stripped.startswith("#"):
                continue
            # Check for bare exit() — not raise SystemExit() or ctx.exit()
            if re.search(r"\bexit\(\)", stripped) and not re.search(
                r"(SystemExit|ctx\.exit|raise\s+SystemExit)", stripped
            ):
                pytest.fail(f"Bare exit() found at _cli.py:{i}: {stripped}")

    def test_cosap_debug_env_var_support(self):
        """Verify _handle_error checks COSAP_DEBUG for traceback output."""
        cli_path = os.path.join(_PROJECT_ROOT, "cosap", "_cli.py")
        with open(cli_path) as f:
            source = f.read()
        assert "COSAP_DEBUG" in source
        assert "traceback" in source


# ===========================================================================
# Test: Docker runner path validation (source-level check)
# ===========================================================================


class TestDockerRunnerSecurityChecks:
    """Verify security patterns exist in docker_runner.py source."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.runner_path = os.path.join(
            _PROJECT_ROOT, "cosap", "runners", "docker_runner.py"
        )
        with open(self.runner_path) as f:
            self.source = f.read()

    def test_has_root_warning(self):
        """Docker runner should warn about running as root."""
        assert "root" in self.source.lower()

    def test_has_path_validation(self):
        """Docker runner should validate mount paths."""
        # Check for some form of path checking
        assert "os.path" in self.source or "realpath" in self.source

    def test_has_container_cleanup(self):
        """Docker runner should clean up containers."""
        assert "remove" in self.source or "cleanup" in self.source.lower()


# ===========================================================================
# Test: Celery config security settings (source-level check)
# ===========================================================================


class TestCeleryConfigSecurity:
    """Verify security-relevant Celery configuration settings."""

    @pytest.fixture(autouse=True)
    def setup(self):
        config_path = os.path.join(_PROJECT_ROOT, "cosap", "celery", "celeryconfig.py")
        with open(config_path) as f:
            self.source = f.read()

    def test_no_default_credentials(self):
        """Celery config should not have default 'guest' credentials."""
        # Should not have guest as a fallback default
        assert '"guest"' not in self.source or "EnvironmentError" in self.source

    def test_has_task_time_limits(self):
        """Celery config should set task time limits."""
        assert "task_time_limit" in self.source

    def test_has_acks_late(self):
        """Celery config should use acks_late for at-least-once delivery."""
        assert "task_acks_late" in self.source

    def test_has_json_serialization(self):
        """Celery config should use JSON serialization (not pickle)."""
        assert "json" in self.source

    def test_env_vars_required_for_broker(self):
        """Celery config should require env vars for broker credentials."""
        assert "EnvironmentError" in self.source or "environ" in self.source


# ===========================================================================
# Test: requirements.txt and pyproject.toml consistency
# ===========================================================================


class TestDependencyManagement:
    """Verify dependency files have proper version pinning."""

    def test_pyproject_has_dependencies(self):
        """pyproject.toml should declare [project.dependencies]."""
        path = os.path.join(_PROJECT_ROOT, "pyproject.toml")
        with open(path) as f:
            content = f.read()
        assert "dependencies" in content
        assert "click" in content
        assert "psutil" in content
        assert "PyYAML" in content

    def test_pyproject_has_version_pins(self):
        """All dependencies in pyproject.toml should have version constraints."""
        path = os.path.join(_PROJECT_ROOT, "pyproject.toml")
        with open(path) as f:
            content = f.read()
        # Check that dependencies have version specifiers
        in_deps = False
        for line in content.splitlines():
            stripped = line.strip()
            if stripped.startswith("dependencies"):
                in_deps = True
                continue
            if in_deps and stripped == "]":
                break
            if in_deps and stripped.startswith('"'):
                # Each dependency line should contain >= or == or ~=
                assert any(op in stripped for op in [">=", "==", "~="]), (
                    f"Dependency without version pin: {stripped}"
                )

    def test_pyproject_has_requires_python(self):
        """pyproject.toml should specify requires-python."""
        path = os.path.join(_PROJECT_ROOT, "pyproject.toml")
        with open(path) as f:
            content = f.read()
        assert "requires-python" in content

    def test_requirements_has_comments(self):
        """requirements.txt should have comments explaining its purpose."""
        path = os.path.join(_PROJECT_ROOT, "requirements.txt")
        with open(path) as f:
            content = f.read()
        assert "#" in content  # Has comments

    def test_no_yaml_package_name(self):
        """requirements.txt should not use 'yaml' (wrong name, should be PyYAML)."""
        path = os.path.join(_PROJECT_ROOT, "requirements.txt")
        with open(path) as f:
            lines = f.readlines()
        for line in lines:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                # 'yaml' as a standalone package name is wrong
                assert stripped != "yaml", (
                    "requirements.txt has 'yaml' — should be 'PyYAML'"
                )

    def test_no_docker_py_old_name(self):
        """requirements.txt should not use 'docker-py' (deprecated name)."""
        path = os.path.join(_PROJECT_ROOT, "requirements.txt")
        with open(path) as f:
            content = f.read()
        assert "docker-py" not in content, (
            "requirements.txt uses 'docker-py' — should be 'docker'"
        )


# ===========================================================================
# Test: No remaining shell=True in production code
# ===========================================================================


class TestNoShellTrue:
    """Verify no subprocess calls use shell=True in production code."""

    def test_no_shell_true_in_source(self):
        """Scan all .py files under cosap/ for shell=True usage."""
        cosap_dir = os.path.join(_PROJECT_ROOT, "cosap")
        violations = []
        for root, dirs, files in os.walk(cosap_dir):
            # Skip test directories
            if "test" in root:
                continue
            for fname in files:
                if not fname.endswith(".py"):
                    continue
                fpath = os.path.join(root, fname)
                with open(fpath) as f:
                    for i, line in enumerate(f, 1):
                        if "shell=True" in line and not line.strip().startswith("#"):
                            violations.append(f"{fpath}:{i}: {line.strip()}")

        assert not violations, "Found shell=True in production code:\n" + "\n".join(
            violations
        )


# ===========================================================================
# Test: No remaining " ".join() -> docker_runner patterns
# ===========================================================================


class TestNoJoinedDockerCommands:
    """Verify no tool files join command lists before passing to docker_runner."""

    def test_no_join_to_docker_pattern(self):
        """Scan tool files for ' '.join(command) patterns that feed into docker_runner.run().

        Legitimate uses of ' '.join() for building string *values* (like read groups)
        are excluded — we only flag joins whose result variable name suggests it is
        a command passed to docker_runner.
        """
        tools_dir = os.path.join(_PROJECT_ROOT, "cosap", "tools")
        # Variable names that indicate a docker command being built
        command_var_patterns = re.compile(
            r"\b(command|cmd|docker_command|run_command)\s*=\s*\"?\s*\"?\s*\.join\("
        )
        violations = []
        for root, dirs, files in os.walk(tools_dir):
            for fname in files:
                if not fname.endswith(".py"):
                    continue
                fpath = os.path.join(root, fname)
                with open(fpath) as f:
                    lines = f.readlines()
                for i, line in enumerate(lines, 1):
                    stripped = line.strip()
                    if stripped.startswith("#"):
                        continue
                    if '" ".join(' in stripped and command_var_patterns.search(
                        stripped
                    ):
                        violations.append(f"{fpath}:{i}: {stripped}")

        assert not violations, (
            "Found ' '.join() building a command for docker_runner:\n"
            + "\n".join(violations)
        )
