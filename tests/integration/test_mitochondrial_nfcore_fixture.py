import gzip
import json
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

import pytest

RUN_ENV = "COSAP_RUN_MITO_NFCORE_TEST"


def _require_integration_environment():
    if os.getenv(RUN_ENV) != "1":
        pytest.skip(f"Set {RUN_ENV}=1 to run the nf-core mitochondrial fixture test.")
    pytest.importorskip("pysam")
    for executable in ("bwa", "gatk", "samtools", "snakemake"):
        if shutil.which(executable) is None:
            pytest.skip(f"{executable} is required for the mitochondrial fixture test.")
    if os.getenv("COSAP_CONTAINER_RUNTIME", "docker") == "docker" and shutil.which("docker") is None:
        pytest.skip("Docker is required when COSAP_CONTAINER_RUNTIME=docker.")


def _manifest() -> dict:
    manifest_path = Path(__file__).parents[1] / "data" / "mitochondrial" / "manifest.json"
    return json.loads(manifest_path.read_text())


def _download_fixture_files(manifest: dict, root: Path) -> None:
    for entry in manifest["files"]:
        destination = root / entry["destination"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and destination.stat().st_size == entry["size"]:
            continue
        urllib.request.urlretrieve(entry["url"], destination)
        actual_size = destination.stat().st_size
        assert actual_size == entry["size"], (
            f"Downloaded size mismatch for {destination}: "
            f"expected {entry['size']}, observed {actual_size}"
        )


def _run(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def _stage_cosap_library(root: Path) -> Path:
    fasta_gz = root / "Homo_sapiens_assembly38.fasta.gz"
    fasta = root / "Homo_sapiens_assembly38.fasta"
    if not fasta.exists():
        with gzip.open(fasta_gz, "rb") as input_file, fasta.open("wb") as output_file:
            shutil.copyfileobj(input_file, output_file)

    if not (root / "Homo_sapiens_assembly38.dict").exists():
        _run(
            [
                "gatk",
                "CreateSequenceDictionary",
                "-R",
                str(fasta),
                "-O",
                str(root / "Homo_sapiens_assembly38.dict"),
            ],
            root,
        )
    if not (root / "Homo_sapiens_assembly38.fasta.bwt").exists():
        _run(["bwa", "index", str(fasta)], root)
    if not (root / "Homo_sapiens_assembly38.fasta.fai").exists():
        _run(["samtools", "faidx", str(fasta)], root)
    return root


def _configure_library(monkeypatch, library_path: Path) -> None:
    from cosap._config import AppConfig
    from cosap._library_paths import _LibraryMeta

    monkeypatch.setenv("COSAP_LIBRARY_PATH", str(library_path))
    monkeypatch.setenv("COSAP_REF_VERSION", "hg38")
    AppConfig.LIBRARY_PATH = str(library_path)
    AppConfig.REF_VERSION = "hg38"
    _LibraryMeta._instances.clear()


def _assert_vcf_contract(vcf_path: Path) -> None:
    import pysam

    assert vcf_path.exists()
    assert vcf_path.stat().st_size > 0
    assert Path(f"{vcf_path}.tbi").exists() or Path(f"{vcf_path}.idx").exists()

    with pysam.VariantFile(str(vcf_path)) as vcf:
        assert "chrM" in vcf.header.contigs
        for record in vcf:
            assert record.contig == "chrM"


@pytest.mark.integration
def test_mitochondrial_pipeline_with_nfcore_fixture(tmp_path, monkeypatch):
    _require_integration_environment()
    manifest = _manifest()
    fixture_root = tmp_path / "nfcore_variantcatalogue_fixture"
    _download_fixture_files(manifest, fixture_root)
    library_path = _stage_cosap_library(fixture_root)
    _configure_library(monkeypatch, library_path)

    from cosap.workflows import MitochondrialPipeline, MitochondrialPipelineInput

    workdir = tmp_path / "workdir"
    workdir.mkdir()
    pipeline = MitochondrialPipeline(
        MitochondrialPipelineInput(
            WORKDIR=str(workdir),
            SAMPLE_NAME="NA12718",
            FASTQS=(
                str(fixture_root / "fastq" / "NA12718_R1.fastq.gz"),
                str(fixture_root / "fastq" / "NA12718_R2.fastq.gz"),
            ),
            SKIP_RECALIBRATION=True,
            RUN_HAPLOCHECKER=False,
        )
    )
    pipeline.run_pipeline()

    output_dir = workdir / "MITOCHONDRIAL_VARIANTS" / "gatk_mitochondrial"
    expected_outputs = [
        "NA12718.chrM.bam",
        "NA12718.chrM.aligned.bam",
        "NA12718.chrM.shifted.aligned.bam",
        "NA12718.chrM.raw.vcf.gz",
        "NA12718.chrM.shifted.raw.vcf.gz",
        "NA12718.chrM.merged.raw.vcf.gz",
        "NA12718.chrM.merged.stats",
        "NA12718.chrM.vcf.gz",
        "NA12718.chrM.split.vcf.gz",
        "NA12718.chrM.coverage.metrics",
        "NA12718.chrM.duplicate.metrics",
    ]
    for relative_path in expected_outputs:
        output_path = output_dir / relative_path
        assert output_path.exists(), f"Missing mitochondrial output: {output_path}"
        assert output_path.stat().st_size > 0, f"Empty mitochondrial output: {output_path}"

    _assert_vcf_contract(output_dir / "NA12718.chrM.vcf.gz")
    _assert_vcf_contract(output_dir / "NA12718.chrM.split.vcf.gz")
