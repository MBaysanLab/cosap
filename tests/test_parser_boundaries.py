import subprocess
import sys

from cosap_formats import AnnotatedCnvVcfParser

from cosap import _utils
from cosap.parsers.annotated_cnv_vcf_parser import (
    AnnotatedCnvVcfParser as CompatibilityAnnotatedCnvVcfParser,
)
from cosap.utils import vcf_utils


def test_vcf_utils_fast_conversion_delegates_to_cosap_formats(monkeypatch):
    expected = '[{"CHROM":"1","POS":42}]'
    calls = []

    def fake_parse(path, caller_type, sample_name, filter_status):
        calls.append((path, caller_type, sample_name, filter_status))
        return expected

    monkeypatch.setattr(vcf_utils, "parse_vcf_json", fake_parse)

    result = vcf_utils.VCFUtils.convert_vcf_to_json_fast(
        "variants.vcf", "mutect2", "tumor", "PASS"
    )

    assert result is expected
    assert calls == [("variants.vcf", "mutect2", "tumor", "PASS")]


def test_legacy_utils_fast_conversion_preserves_delegate_output(monkeypatch):
    expected = '[{"CHROM":"2","POS":84}]'
    monkeypatch.setattr(
        vcf_utils.VCFUtils,
        "convert_vcf_to_json_fast",
        lambda *args: expected,
    )

    assert (
        _utils.convert_vcf_to_json_fast("variants.vcf", sample_name="normal")
        is expected
    )


def test_annotated_cnv_parser_is_compatibility_reexport():
    assert CompatibilityAnnotatedCnvVcfParser is AnnotatedCnvVcfParser


def test_parser_compatibility_imports_do_not_require_celery():
    code = """
import builtins

real_import = builtins.__import__

def import_without_celery(name, *args, **kwargs):
    if name == "celery" or name.startswith("celery."):
        raise AssertionError("parser compatibility imports must not load Celery")
    return real_import(name, *args, **kwargs)

builtins.__import__ = import_without_celery

from cosap._utils import convert_vcf_to_json_fast
from cosap.parsers.annotated_cnv_vcf_parser import AnnotatedCnvVcfParser
from cosap_formats import AnnotatedCnvVcfParser as FormatsAnnotatedCnvVcfParser

assert AnnotatedCnvVcfParser is FormatsAnnotatedCnvVcfParser
assert callable(convert_vcf_to_json_fast)
"""

    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
