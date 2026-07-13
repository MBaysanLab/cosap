from unittest.mock import Mock

import pytest

pytest.importorskip("psutil")

from cosap import _config
from cosap._pipeline_config import MappingKeys
from cosap.tools.mappers._bowtie_mapper import Bowtie2Mapper
from cosap.tools.mappers._bwa2_mapper import BWA2Mapper
from cosap.tools.mappers._bwa_mapper import BWAMapper
from cosap.tools.mappers._mappers import _Mapper


@pytest.fixture
def mapper_config(tmp_path):
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    r1.write_text("@INST:1:FC:1:1:1:1:AAAA 1:N:0:ATCG\n")
    r2.write_text("@INST:1:FC:1:1:1:1:AAAA 2:N:0:ATCG\n")

    return {
        MappingKeys.INPUT: {"1": str(r1), "2": str(r2)},
        MappingKeys.OUTPUT: str(output_dir / "out.bam"),
        MappingKeys.OUTPUT_DIR: str(output_dir),
        MappingKeys.PARAMS: {MappingKeys.READ_GROUP: {MappingKeys.RG_SM: "sample"}},
        MappingKeys.POST_PROCESSING: False,
        MappingKeys.WORKDIR: str(output_dir),
    }


def test_samtools_sort_command_uses_configured_threads_and_memory(monkeypatch):
    app_config = _config.AppConfig()
    monkeypatch.setattr(app_config, "SAMTOOLS_SORT_THREADS", 2)
    monkeypatch.setattr(app_config, "SAMTOOLS_SORT_MEM_PER_THREAD", "512M")

    command = _Mapper._samtools_sort_command(app_config, "out.bam")

    assert command[:6] == ["samtools", "sort", "-@", "2", "-m", "512M"]


def test_bwa_parabricks_command_splits_read_group_arguments(mapper_config):
    mapper_config[MappingKeys.OUTPUT] = "out.bam"
    mapper_config[MappingKeys.OUTPUT_DIR] = "output"
    mapper_config[MappingKeys.PARAMS][MappingKeys.READ_GROUP] = {
        MappingKeys.RG_ID: 0,
        MappingKeys.RG_SM: "normal",
        MappingKeys.RG_LB: 0,
        MappingKeys.RG_PL: "ILLUMINA",
    }

    read_group = BWAMapper._create_read_group(
        mapper_config=mapper_config, is_parabricks=True
    )
    command = BWAMapper._create_parabricks_command(
        mapper_config=mapper_config,
        read_group=read_group,
        library_paths=Mock(BWA_ASSEMBLY="ref.fasta"),
    )

    assert command[command.index("--out-bam") + 1] == "output/out.bam"
    assert "--bwa-options=-K 10000000" in command
    assert "--bwa-options=-'K 10000000'" not in command
    assert "--read-group-id-prefix 0 --read-group-sm normal" not in command
    assert command[-8:] == [
        "--read-group-id-prefix",
        "0",
        "--read-group-sm",
        "normal",
        "--read-group-lb",
        "0",
        "--read-group-pl",
        "ILLUMINA",
    ]


@pytest.mark.parametrize("mapper", [BWAMapper, BWA2Mapper, Bowtie2Mapper])
def test_mapper_streams_with_popen(monkeypatch, mapper, mapper_config):
    monkeypatch.setattr(_config.AppConfig, "MAX_THREADS_PER_JOB", 2)
    monkeypatch.setattr(_config.AppConfig, "SAMTOOLS_SORT_THREADS", 2)
    monkeypatch.setattr(_config.AppConfig, "SAMTOOLS_SORT_MEM_PER_THREAD", "512M")

    class _Stdout:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    first = Mock()
    first.stdout = _Stdout()
    first.communicate.return_value = (b"", b"")
    first.returncode = 0

    second = Mock()
    second.communicate.return_value = (b"", b"")
    second.returncode = 0

    popen_calls = []

    def _fake_popen(*args, **kwargs):
        popen_calls.append((args, kwargs))
        return first if len(popen_calls) == 1 else second

    module = __import__(mapper.__module__, fromlist=["Popen"])
    monkeypatch.setattr(module, "Popen", _fake_popen)

    mapper.map(mapper_config)

    assert len(popen_calls) == 2
    assert popen_calls[0][1]["stdout"] is not None
    assert popen_calls[1][1]["stdin"] is first.stdout
    assert first.stdout.closed is True


@pytest.mark.parametrize("mapper", [BWAMapper, BWA2Mapper, Bowtie2Mapper])
def test_mapper_raises_when_sort_fails(monkeypatch, mapper, mapper_config):
    monkeypatch.setattr(_config.AppConfig, "MAX_THREADS_PER_JOB", 2)
    monkeypatch.setattr(_config.AppConfig, "SAMTOOLS_SORT_THREADS", 2)
    monkeypatch.setattr(_config.AppConfig, "SAMTOOLS_SORT_MEM_PER_THREAD", "512M")

    class _Stdout:
        def close(self):
            return None

    first = Mock()
    first.stdout = _Stdout()
    first.communicate.return_value = (b"", b"")
    first.returncode = 0

    second = Mock()
    second.communicate.return_value = (b"", b"sort failed")
    second.returncode = 1

    popen_calls = []

    def _fake_popen(*args, **kwargs):
        popen_calls.append((args, kwargs))
        return first if len(popen_calls) == 1 else second

    module = __import__(mapper.__module__, fromlist=["Popen"])
    monkeypatch.setattr(module, "Popen", _fake_popen)

    with pytest.raises(Exception, match="samtools sort command"):
        mapper.map(mapper_config)
