from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from cosap.tools.annotators._exomiser import Exomiser


@pytest.mark.parametrize("trailing_separator", [False, True])
def test_create_job_file_joins_output_directory(
    monkeypatch, tmp_path, trailing_separator
):
    variant_file = SimpleNamespace(
        header=SimpleNamespace(samples=["sample_1"])
    )
    monkeypatch.setattr(
        "cosap.tools.annotators._exomiser.pysam.VariantFile",
        lambda _path: variant_file,
    )
    output_dir = str(tmp_path) + ("/" if trailing_separator else "")

    job_file = Exomiser._create_job_file(
        "input.vcf.gz",
        ["HP:0000001"],
        str(tmp_path / "annotated_sample_1"),
        output_dir,
    )

    assert Path(job_file) == tmp_path / "exomiser_job_sample_1.yaml"
    assert Path(job_file).is_file()
    assert yaml.safe_load(Path(job_file).read_text())["phenopacket"]["subject"] == {
        "id": "sample_1"
    }
