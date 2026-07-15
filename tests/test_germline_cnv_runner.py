from pathlib import Path

from cosap.tools.cnv_callers._gatk_germline_cnv_caller import GATKGermlineCNVCaller


def _write_interval_list(path: Path, sizes):
    with path.open("w") as handle:
        handle.write("@HD\tVN:1.6\n")
        for i, size in enumerate(sizes, start=1):
            handle.write(f"chr1\t{i * 1000}\t{i * 1000 + size - 1}\t+\ttarget_{i}\n")


def test_cram_sample_name_strips_cram_suffix():
    assert (
        GATKGermlineCNVCaller._sample_name_from_alignment_path("/data/HG00512.cram")
        == "HG00512"
    )
    assert (
        GATKGermlineCNVCaller._sample_name_from_alignment_path("/data/HG00512.bam")
        == "HG00512"
    )


def test_resolve_intervals_per_shard_targets_base_pairs_for_short_exome_targets(tmp_path):
    interval_list = tmp_path / "targets.interval_list"
    _write_interval_list(interval_list, [200] * 100_000)

    intervals_per_shard = GATKGermlineCNVCaller._resolve_intervals_per_shard(
        str(interval_list), {"target_bases_per_shard": 10_000_000}
    )

    assert intervals_per_shard == 50_000


def test_resolve_intervals_per_shard_keeps_explicit_override(tmp_path):
    interval_list = tmp_path / "targets.interval_list"
    _write_interval_list(interval_list, [200] * 100_000)

    intervals_per_shard = GATKGermlineCNVCaller._resolve_intervals_per_shard(
        str(interval_list), {"intervals_per_shard": 5000}
    )

    assert intervals_per_shard == 5000
