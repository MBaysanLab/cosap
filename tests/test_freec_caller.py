from pathlib import Path

import pytest

from cosap._docker_images import DockerImages
from cosap._pipeline_config import CNVCallingKeys
from cosap.tools.cnv_callers import FREEC, CNVCallerFactory


def _touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
    return str(path)


def _base_config(tmp_path: Path) -> dict:
    bam = _touch(tmp_path / "sample.bam")
    bed = _touch(tmp_path / "targets.bed")
    fai = _touch(tmp_path / "ref.fa.fai")
    chr_files = tmp_path / "chr_files"
    chr_files.mkdir()
    return {
        CNVCallingKeys.TUMOR_INPUT: bam,
        CNVCallingKeys.OUTPUT_DIR: str(tmp_path / "out"),
        CNVCallingKeys.OUTPUT: "segments.tsv",
        CNVCallingKeys.PARAMS: {
            "capture_regions": bed,
            "chr_len_file": fai,
            "chr_files": str(chr_files),
            "ploidy": "2",
            "threads": 2,
        },
    }


def test_factory_creates_freec_case_insensitively():
    assert CNVCallerFactory.create("FREEC") is FREEC


def test_config_for_control_free_wes_includes_chr_files(tmp_path):
    config = _base_config(tmp_path)
    output_dir = str(tmp_path / "out")
    resolved = FREEC._validate_inputs(config)

    lines = FREEC._config_lines(resolved, output_dir)
    text = "\n".join(lines)

    assert "window = 0" in text
    assert "noisyData = TRUE" in text
    assert "readCountThreshold = 50" in text
    assert f"chrFiles = {resolved['chr_files']}" in text
    assert f"captureRegions = {resolved['capture_regions']}" in text
    assert "[control]" not in text


def test_config_for_paired_wes_includes_control(tmp_path):
    config = _base_config(tmp_path)
    config[CNVCallingKeys.NORMAL_INPUT] = _touch(tmp_path / "normal.bam")
    resolved = FREEC._validate_inputs(config)

    text = "\n".join(FREEC._config_lines(resolved, str(tmp_path / "out")))

    assert "[control]" in text
    assert f"mateFile = {resolved['control_input']}" in text


def test_control_free_wes_requires_gc_reference(tmp_path):
    config = _base_config(tmp_path)
    config[CNVCallingKeys.PARAMS].pop("chr_files")

    with pytest.raises(ValueError, match="control-free WES requires"):
        FREEC._validate_inputs(config)


def test_freec_rejects_cram_input(tmp_path):
    config = _base_config(tmp_path)
    config[CNVCallingKeys.TUMOR_INPUT] = _touch(tmp_path / "sample.cram")

    with pytest.raises(ValueError, match="CRAM"):
        FREEC._validate_inputs(config)


def test_normalise_cnvs_writes_stable_header(tmp_path):
    raw = tmp_path / "sample_CNVs"
    raw.write_text("chr1\t10\t20\t1\tloss\nchr2 30 40 3 gain extra\n")
    output = tmp_path / "segments.tsv"

    FREEC._normalise_cnvs(str(raw), str(output))

    assert output.read_text().splitlines() == [
        "chrom\tstart\tend\tcopy_number\tstatus\textra",
        "chr1\t10\t20\t1\tloss\t",
        "chr2\t30\t40\t3\tgain\textra",
    ]


def test_parse_info_accepts_equals_and_tab_delimiters(tmp_path):
    info = tmp_path / "sample_info.txt"
    info.write_text("Good_Polynomial_Fit=True\nCGcontent_Used\tFalse\n")

    assert FREEC._parse_info(str(info)) == {
        "Good_Polynomial_Fit": "True",
        "CGcontent_Used": "False",
    }


def test_control_free_wes_rejects_output_without_gc_normalization(tmp_path):
    config = _base_config(tmp_path)
    out_dir = Path(config[CNVCallingKeys.OUTPUT_DIR])
    out_dir.mkdir(parents=True)
    (out_dir / "sample_info.txt").write_text(
        "Good_Polynomial_Fit\tTrue\nCGcontent_Used\tFalse\n"
    )
    (out_dir / "sample_CNVs").write_text("chr1\t10\t20\t1\tloss\n")

    with pytest.raises(RuntimeError, match="CGcontent_Used=False"):
        FREEC._validate_run_outputs(config, str(out_dir))


def test_paired_wes_allows_output_without_gc_normalization(tmp_path):
    config = _base_config(tmp_path)
    config[CNVCallingKeys.NORMAL_INPUT] = _touch(tmp_path / "normal.bam")
    out_dir = Path(config[CNVCallingKeys.OUTPUT_DIR])
    out_dir.mkdir(parents=True)
    (out_dir / "sample_info.txt").write_text(
        "Good_Polynomial_Fit\tTrue\nCGcontent_Used\tFalse\n"
    )
    raw_cnvs = out_dir / "sample_CNVs"
    raw_cnvs.write_text("chr1\t10\t20\t1\tloss\n")

    assert FREEC._validate_run_outputs(config, str(out_dir)) == str(raw_cnvs)


def test_call_runs_container_and_normalises_output(monkeypatch, tmp_path):
    config = _base_config(tmp_path)
    out_dir = Path(config[CNVCallingKeys.OUTPUT_DIR])

    class FakeRunner:
        def __init__(self, device="cpu"):
            self.device = device

        def run(self, image, command, workdir=None, paths_to_bind=None):
            assert image == DockerImages.FREEC
            assert command[0:2] == ["freec", "-conf"]
            assert workdir == str(out_dir.resolve())
            assert paths_to_bind
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "sample_info.txt").write_text("Good_Polynomial_Fit=True\n")
            (out_dir / "sample_CNVs").write_text("chr1\t10\t20\t1\tloss\n")
            return 0

    monkeypatch.setattr(
        "cosap.tools.cnv_callers._freec_cnv_caller.DockerRunner", FakeRunner
    )

    result = FREEC.call(config)

    assert result == str(out_dir.resolve() / "segments.tsv")
    assert Path(result).read_text().splitlines()[1] == "chr1\t10\t20\t1\tloss\t"
