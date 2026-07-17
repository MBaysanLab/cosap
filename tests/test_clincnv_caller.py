from pathlib import Path

import pytest

from cosap._docker_images import DockerImages
from cosap._pipeline_config import CNVCallingKeys
from cosap.tools.cnv_callers import CNVCallerFactory, ClinCNV
from cosap.tools.cnv_callers import _clincnv_cnv_caller as clincnv_module


def _config(tmp_path: Path) -> dict:
    matrix = tmp_path / "coverage.cov"
    matrix.write_text("chr\tstart\tend\tHG00512\tother\nchr1\t1\t10\t20\t21\n")
    bed = tmp_path / "targets.gc.bed"
    bed.write_text("chr1\t1\t10\t0.5\n")
    scripts = tmp_path / "ClinCNV"
    scripts.mkdir()
    (scripts / "clinCNV.R").write_text("")
    return {
        CNVCallingKeys.OUTPUT_DIR: str(tmp_path / "out"),
        CNVCallingKeys.OUTPUT: "segments.tsv",
        CNVCallingKeys.PARAMS: {
            "coverage_matrix": str(matrix),
            "annotated_bed": str(bed),
            "script_dir": str(scripts),
            "sample_name": "HG00512",
            "hg38": "TRUE",
        },
    }


def test_factory_creates_clincnv_case_insensitively():
    assert CNVCallerFactory.create("ClinCNV") is ClinCNV


def test_clincnv_command_targets_one_cohort_sample(tmp_path):
    resolved = ClinCNV._validate_inputs(_config(tmp_path))
    command = ClinCNV._command(resolved, str(tmp_path / "out"))

    assert "--normalSample" in command
    assert command[command.index("--normalSample") + 1] == "HG00512"
    assert command[command.index("--hg38") + 1] == "TRUE"


def test_clincnv_defaults_to_container_script_dir(tmp_path):
    config = _config(tmp_path)
    del config[CNVCallingKeys.PARAMS]["script_dir"]

    resolved = ClinCNV._validate_inputs(config)

    assert resolved["script_dir"] == "/opt/ClinCNV"
    assert resolved["host_script_dir"] is None


def test_clincnv_requires_sample_in_coverage_matrix(tmp_path):
    config = _config(tmp_path)
    config[CNVCallingKeys.PARAMS]["sample_name"] = "missing"

    with pytest.raises(ValueError, match="coverage-matrix column"):
        ClinCNV._validate_inputs(config)


def test_normalise_clincnv_cnvs(tmp_path):
    raw = tmp_path / "sample_cnvs.tsv"
    raw.write_text(
        "##ANALYSISTYPE=CLINCNV_GERMLINE_SINGLE\n"
        "#chr\tstart\tend\tCN_change\tloglikelihood\n"
        "chr1\t10\t20\t1\t60\nchr2\t30\t40\t3\t70\n"
    )
    output = tmp_path / "segments.tsv"

    ClinCNV._normalise_cnvs(str(raw), str(output))

    assert output.read_text().splitlines() == [
        "chrom\tstart\tend\tcopy_number\tstatus\textra",
        "chr1\t10\t20\t1\tloss\tloglikelihood=60",
        "chr2\t30\t40\t3\tgain\tloglikelihood=70",
    ]


def test_clincnv_call_runs_in_docker(monkeypatch, tmp_path):
    config = _config(tmp_path)
    del config[CNVCallingKeys.PARAMS]["script_dir"]
    calls = []

    class FakeDockerRunner:
        def __init__(self, device="cpu"):
            self.device = device

        def run(self, image, command, workdir=None, paths_to_bind=None):
            calls.append(
                {
                    "image": image,
                    "command": command,
                    "workdir": workdir,
                    "paths_to_bind": paths_to_bind,
                }
            )
            raw_dir = Path(workdir) / "normal" / "HG00512"
            raw_dir.mkdir(parents=True)
            (raw_dir / "HG00512_cnvs.tsv").write_text(
                "#chr\tstart\tend\tCN_change\tloglikelihood\n"
                "chr1\t10\t20\t1\t60\n"
            )
            return 0

    monkeypatch.setattr(clincnv_module, "DockerRunner", FakeDockerRunner)

    output = ClinCNV.call(config)

    assert output == str(tmp_path / "out" / "segments.tsv")
    assert calls[0]["image"] == DockerImages.CLINCNV
    assert calls[0]["command"][0] == "Rscript"
    assert calls[0]["command"][1] == "/opt/ClinCNV/clinCNV.R"
    assert Path(output).read_text().splitlines()[1] == "chr1\t10\t20\t1\tloss\tloglikelihood=60"
