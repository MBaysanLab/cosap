import json
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

import pytest


RUN_ENV = "COSAP_RUN_MITO_NA12878_BENCHMARK"
NA12878_BAM = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/"
    "NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
)
NA12878_1000G_CHRMT_VCF = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
)


def _require_benchmark_environment():
    if os.getenv(RUN_ENV) != "1":
        pytest.skip(f"Set {RUN_ENV}=1 to run the NA12878 1000G mitochondrial benchmark.")
    pytest.importorskip("pysam")
    for executable in ("bcftools", "bwa", "gatk", "samtools", "snakemake"):
        if shutil.which(executable) is None:
            pytest.skip(f"{executable} is required for the mitochondrial benchmark.")
    if os.getenv("COSAP_CONTAINER_RUNTIME", "docker") == "docker" and shutil.which("docker") is None:
        pytest.skip("Docker is required when COSAP_CONTAINER_RUNTIME=docker.")


def _run(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def _run_shell(command: str, cwd: Path) -> None:
    subprocess.run(["bash", "-lc", command], cwd=cwd, check=True)


def _manifest() -> dict:
    manifest_path = Path(__file__).parents[1] / "data" / "mitochondrial" / "manifest.json"
    return json.loads(manifest_path.read_text())


def _download_mitochondrial_references(root: Path) -> Path:
    manifest = _manifest()
    for entry in manifest["files"]:
        if not entry["destination"].startswith("chrM/"):
            continue
        destination = root / entry["destination"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and destination.stat().st_size == entry["size"]:
            continue
        urllib.request.urlretrieve(entry["url"], destination)
        assert destination.stat().st_size == entry["size"]
    return root


def _configure_library(monkeypatch, library_path: Path) -> None:
    from cosap._config import AppConfig
    from cosap._library_paths import _LibraryMeta, _LibraryPaths38

    monkeypatch.setenv("COSAP_LIBRARY_PATH", str(library_path))
    monkeypatch.setenv("COSAP_REF_VERSION", "hg38")
    AppConfig.LIBRARY_PATH = str(library_path)
    AppConfig.REF_VERSION = "hg38"
    mt_dir = library_path / "chrM"
    mitochondrial_paths = {
        "REF_DIR": library_path,
        "MT_DIR": mt_dir,
        "MT_FASTA": mt_dir / "Homo_sapiens_assembly38.chrM.fasta",
        "MT_DICT": mt_dir / "Homo_sapiens_assembly38.chrM.dict",
        "MT_FASTA_INDEX": mt_dir / "Homo_sapiens_assembly38.chrM.fasta.fai",
        "MT_FASTA_AMB": mt_dir / "Homo_sapiens_assembly38.chrM.fasta.amb",
        "MT_FASTA_ANN": mt_dir / "Homo_sapiens_assembly38.chrM.fasta.ann",
        "MT_FASTA_BWT": mt_dir / "Homo_sapiens_assembly38.chrM.fasta.bwt",
        "MT_FASTA_PAC": mt_dir / "Homo_sapiens_assembly38.chrM.fasta.pac",
        "MT_FASTA_SA": mt_dir / "Homo_sapiens_assembly38.chrM.fasta.sa",
        "MT_SHIFTED_FASTA": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta",
        "MT_SHIFTED_DICT": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict",
        "MT_SHIFTED_FASTA_INDEX": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai",
        "MT_SHIFTED_FASTA_AMB": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb",
        "MT_SHIFTED_FASTA_ANN": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann",
        "MT_SHIFTED_FASTA_BWT": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt",
        "MT_SHIFTED_FASTA_PAC": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac",
        "MT_SHIFTED_FASTA_SA": mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa",
        "MT_SHIFT_BACK_CHAIN": mt_dir / "ShiftBack.chain",
        "MT_BLACKLIST": mt_dir / "blacklist_sites.hg38.chrM.bed",
        "MT_BLACKLIST_INDEX": mt_dir / "blacklist_sites.hg38.chrM.bed.idx",
    }
    for key, value in mitochondrial_paths.items():
        setattr(_LibraryPaths38, key, str(value))
        _LibraryPaths38.__dataclass_fields__[key].default = str(value)
    _LibraryMeta._instances.clear()


def _stage_na12878_chrM_bam(workdir: Path) -> tuple[Path, Path]:
    raw_bam = workdir / "NA12878.MT.1000g.bam"
    chrM_bam = workdir / "NA12878.chrM.1000g.bam"
    header = workdir / "NA12878.chrM.header.sam"

    if not raw_bam.exists():
        _run(["samtools", "view", "-bh", NA12878_BAM, "MT", "-o", str(raw_bam)], workdir)

    if not chrM_bam.exists():
        original_header = subprocess.check_output(
            ["samtools", "view", "-H", str(raw_bam)], cwd=workdir, text=True
        )
        header.write_text(original_header.replace("\tSN:MT\t", "\tSN:chrM\t"))
        _run_shell(
            f"samtools reheader {header} {raw_bam} > {chrM_bam}",
            workdir,
        )

    if not Path(f"{chrM_bam}.bai").exists():
        _run(["samtools", "index", str(chrM_bam)], workdir)
    return chrM_bam, Path(f"{chrM_bam}.bai")


def _normalise_1000g_truth(workdir: Path, reference: Path) -> Path:
    selected = workdir / "NA12878.1000g.chrMT.selected.vcf.gz"
    renamed = workdir / "NA12878.1000g.chrM.selected.vcf.gz"
    normalised = workdir / "NA12878.1000g.chrM.normalised.vcf.gz"
    rename_map = workdir / "rename_chromosomes.tsv"
    rename_map.write_text("MT\tchrM\n")

    if not selected.exists():
        _run(
            [
                "bcftools",
                "view",
                "-s",
                "NA12878",
                "-i",
                'GT!="0" && GT!="0/0" && GT!="./."',
                "-Oz",
                "-o",
                str(selected),
                NA12878_1000G_CHRMT_VCF,
            ],
            workdir,
        )
        _run(["bcftools", "index", "-t", str(selected)], workdir)

    if not normalised.exists():
        _run_shell(
            "bcftools annotate "
            f"--rename-chrs {rename_map} {selected} -Ou | "
            f"bcftools norm -m-any -f {reference} -Oz -o {renamed}",
            workdir,
        )
        _run(["bcftools", "index", "-t", str(renamed)], workdir)
        _run(
            [
                "bcftools",
                "view",
                "-i",
                'GT!="0" && GT!="0/0" && GT!="./."',
                "-Oz",
                "-o",
                str(normalised),
                str(renamed),
            ],
            workdir,
        )
        _run(["bcftools", "index", "-t", str(normalised)], workdir)
    return normalised


def _mitochondrial_reference_files(root: Path) -> dict[str, str]:
    mt_dir = root / "chrM"
    return {
        "standard FASTA": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta"),
        "standard FASTA index": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta.fai"),
        "standard dictionary": str(mt_dir / "Homo_sapiens_assembly38.chrM.dict"),
        "standard BWA .amb": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta.amb"),
        "standard BWA .ann": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta.ann"),
        "standard BWA .bwt": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta.bwt"),
        "standard BWA .pac": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta.pac"),
        "standard BWA .sa": str(mt_dir / "Homo_sapiens_assembly38.chrM.fasta.sa"),
        "shifted FASTA": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"),
        "shifted FASTA index": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai"),
        "shifted dictionary": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict"),
        "shifted BWA .amb": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb"),
        "shifted BWA .ann": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann"),
        "shifted BWA .bwt": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt"),
        "shifted BWA .pac": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac"),
        "shifted BWA .sa": str(mt_dir / "Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa"),
        "shift-back chain": str(mt_dir / "ShiftBack.chain"),
        "blacklist": str(mt_dir / "blacklist_sites.hg38.chrM.bed"),
        "blacklist index": str(mt_dir / "blacklist_sites.hg38.chrM.bed.idx"),
    }


def _normalise_calls(workdir: Path, reference: Path, calls: Path) -> Path:
    normalised = workdir / "NA12878.cosap.chrM.normalised.vcf.gz"
    if not normalised.exists():
        _run(
            [
                "bcftools",
                "norm",
                "-m-any",
                "-f",
                str(reference),
                "-Oz",
                "-o",
                str(normalised),
                str(calls),
            ],
            workdir,
        )
        _run(["bcftools", "index", "-t", str(normalised)], workdir)
    return normalised


def _variant_keys(vcf_path: Path, pass_only: bool) -> set[tuple[str, int, str, str]]:
    import pysam

    keys = set()
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            if pass_only and set(record.filter.keys()) not in (set(), {"PASS"}):
                continue
            for alt in record.alts or []:
                keys.add((record.contig, record.pos, record.ref, alt))
    return keys


def _write_metrics(workdir: Path, truth_vcf: Path, calls_vcf: Path) -> dict:
    truth = _variant_keys(truth_vcf, pass_only=False)
    calls = _variant_keys(calls_vcf, pass_only=True)
    true_positives = truth & calls
    false_positives = calls - truth
    false_negatives = truth - calls
    precision = len(true_positives) / len(calls) if calls else 0.0
    recall = len(true_positives) / len(truth) if truth else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    metrics = {
        "benchmark": "COSAP NA12878 mitochondrial PASS calls vs 1000 Genomes Phase 3 callMom chrMT calls",
        "truth_vcf": str(truth_vcf),
        "calls_vcf": str(calls_vcf),
        "truth_count": len(truth),
        "pass_call_count": len(calls),
        "true_positive_count": len(true_positives),
        "false_positive_count": len(false_positives),
        "false_negative_count": len(false_negatives),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "true_positives": sorted(true_positives),
        "false_positives": sorted(false_positives),
        "false_negatives": sorted(false_negatives),
    }
    (workdir / "benchmark_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))
    return metrics


@pytest.mark.integration
@pytest.mark.benchmark
def test_na12878_1000g_mitochondrial_concordance(tmp_path, monkeypatch):
    _require_benchmark_environment()
    fixture_root = _download_mitochondrial_references(tmp_path / "mitochondrial_references")
    _configure_library(monkeypatch, fixture_root)
    input_bam, input_bai = _stage_na12878_chrM_bam(tmp_path)

    from cosap._formats import FileFormats
    from cosap._pipeline_config import MitochondrialCallingKeys
    from cosap.tools.mitochondrial_callers import GATKMitochondrialCaller

    monkeypatch.setattr(
        GATKMitochondrialCaller,
        "_reference_files",
        classmethod(lambda cls, paths: _mitochondrial_reference_files(fixture_root)),
    )

    output_dir = tmp_path / "MITOCHONDRIAL_VARIANTS" / "gatk_mitochondrial"
    output_dir.mkdir(parents=True)
    name = "NA12878"
    GATKMitochondrialCaller.call(
        {
            MitochondrialCallingKeys.INPUT_BAM: str(input_bam),
            MitochondrialCallingKeys.INPUT_BAM_INDEX: str(input_bai),
            MitochondrialCallingKeys.OUTPUT_DIR: str(output_dir),
            MitochondrialCallingKeys.CONTIG_NAME: "chrM",
            MitochondrialCallingKeys.RUN_HAPLOCHECKER: False,
            MitochondrialCallingKeys.PARAMS: {},
            MitochondrialCallingKeys.SUBSET_BAM: FileFormats.MT_SUBSET_BAM.format(identification=name),
            MitochondrialCallingKeys.UNMAPPED_BAM: FileFormats.MT_UNMAPPED_BAM.format(identification=name),
            MitochondrialCallingKeys.STANDARD_ALIGNED_BAM: FileFormats.MT_STANDARD_ALIGNED_BAM.format(identification=name),
            MitochondrialCallingKeys.SHIFTED_ALIGNED_BAM: FileFormats.MT_SHIFTED_ALIGNED_BAM.format(identification=name),
            MitochondrialCallingKeys.STANDARD_RAW_VCF: FileFormats.MT_STANDARD_RAW_VCF.format(identification=name),
            MitochondrialCallingKeys.SHIFTED_RAW_VCF: FileFormats.MT_SHIFTED_RAW_VCF.format(identification=name),
            MitochondrialCallingKeys.MERGED_RAW_VCF: FileFormats.MT_MERGED_RAW_VCF.format(identification=name),
            MitochondrialCallingKeys.MERGED_STATS: FileFormats.MT_MERGED_STATS.format(identification=name),
            MitochondrialCallingKeys.FILTERED_VCF: FileFormats.MT_FILTERED_VCF.format(identification=name),
            MitochondrialCallingKeys.SPLIT_VCF: FileFormats.MT_SPLIT_VCF.format(identification=name),
            MitochondrialCallingKeys.COVERAGE_METRICS: FileFormats.MT_COVERAGE_METRICS.format(identification=name),
            MitochondrialCallingKeys.DUPLICATE_METRICS: FileFormats.MT_DUPLICATE_METRICS.format(identification=name),
            MitochondrialCallingKeys.CONTAMINATION_REPORT: FileFormats.MT_CONTAMINATION_REPORT.format(identification=name),
        }
    )

    calls = output_dir / "NA12878.chrM.vcf.gz"
    assert calls.exists() and calls.stat().st_size > 0

    reference = fixture_root / "chrM" / "Homo_sapiens_assembly38.chrM.fasta"
    truth_vcf = _normalise_1000g_truth(tmp_path, reference)
    calls_vcf = _normalise_calls(tmp_path, reference, calls)
    metrics = _write_metrics(tmp_path, truth_vcf, calls_vcf)

    assert metrics["truth_count"] > 0
    assert metrics["pass_call_count"] > 0
