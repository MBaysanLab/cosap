"""GATK mitochondrial short-variant discovery workflow."""

import os
import shlex
from pathlib import Path
from typing import Dict, Iterable

from ..._config import AppConfig
from ..._docker_images import DockerImages
from ..._library_paths import LibraryPaths
from ..._pipeline_config import MitochondrialCallingKeys
from ...runners.docker_runner import DockerRunner
from ._mitochondrial_callers import _MitochondrialCaller


class GATKMitochondrialCaller(_MitochondrialCaller):
    """Run GATK's dual-reference mitochondrial Mutect2 workflow on hg38."""

    @staticmethod
    def _path(value: str, output_dir: str) -> str:
        return value if os.path.isabs(value) else os.path.join(output_dir, value)

    @classmethod
    def _reference_files(cls, paths: LibraryPaths) -> Dict[str, str]:
        files = {
            "standard FASTA": paths.MT_FASTA,
            "standard FASTA index": paths.MT_FASTA_INDEX,
            "standard dictionary": paths.MT_DICT,
            "standard BWA .amb": paths.MT_FASTA_AMB,
            "standard BWA .ann": paths.MT_FASTA_ANN,
            "standard BWA .bwt": paths.MT_FASTA_BWT,
            "standard BWA .pac": paths.MT_FASTA_PAC,
            "standard BWA .sa": paths.MT_FASTA_SA,
            "shifted FASTA": paths.MT_SHIFTED_FASTA,
            "shifted FASTA index": paths.MT_SHIFTED_FASTA_INDEX,
            "shifted dictionary": paths.MT_SHIFTED_DICT,
            "shifted BWA .amb": paths.MT_SHIFTED_FASTA_AMB,
            "shifted BWA .ann": paths.MT_SHIFTED_FASTA_ANN,
            "shifted BWA .bwt": paths.MT_SHIFTED_FASTA_BWT,
            "shifted BWA .pac": paths.MT_SHIFTED_FASTA_PAC,
            "shifted BWA .sa": paths.MT_SHIFTED_FASTA_SA,
            "shift-back chain": paths.MT_SHIFT_BACK_CHAIN,
            "blacklist": paths.MT_BLACKLIST,
            "blacklist index": paths.MT_BLACKLIST_INDEX,
        }
        missing = [f"{label}: {path}" for label, path in files.items() if not os.path.isfile(path)]
        if missing:
            raise FileNotFoundError(
                "Missing hg38 mitochondrial reference assets:\n" + "\n".join(missing)
            )
        return files

    @staticmethod
    def _quote(command: Iterable[object]) -> str:
        return " ".join(shlex.quote(str(part)) for part in command)

    @classmethod
    def _run(cls, runner: DockerRunner, command: str, output_dir: str, bind_paths: list) -> None:
        return_code = runner.run(
            image=DockerImages.COSAP,
            command=["bash", "-lc", f"set -euo pipefail; {command}"],
            workdir=output_dir,
            paths_to_bind=bind_paths,
        )
        if return_code not in (0, None):
            raise RuntimeError(f"Mitochondrial workflow command failed with return code {return_code}")

    @classmethod
    def _run_haplochecker(cls, output_dir: str, input_vcf: str, report: str, bind_paths: list) -> None:
        runner = DockerRunner(device="cpu")
        command = (
            "set -euo pipefail; "
            f"mkdir -p haplochecker_input; cp {shlex.quote(input_vcf)} haplochecker_input/input.vcf.gz; "
            "java -jar /usr/mtdnaserver/haplocheckCLI.jar haplochecker_input; "
            f"test -f output && sed 's/\"//g' output > {shlex.quote(report)}"
        )
        return_code = runner.run(
            image=DockerImages.HAPLOCHECKER,
            command=["bash", "-lc", command],
            workdir=output_dir,
            paths_to_bind=bind_paths,
        )
        if return_code not in (0, None):
            raise RuntimeError(f"Haplochecker failed with return code {return_code}")

    @classmethod
    def call(cls, caller_config: dict) -> None:
        if AppConfig.REF_VERSION.lower() != "hg38":
            raise ValueError("GATK mitochondrial calling currently supports hg38 only.")
        paths = LibraryPaths()
        refs = cls._reference_files(paths)
        output_dir = cls._path(
            caller_config[MitochondrialCallingKeys.OUTPUT_DIR],
            caller_config.get(MitochondrialCallingKeys.WORKDIR, os.getcwd()),
        )
        os.makedirs(output_dir, exist_ok=True)

        input_bam = os.path.abspath(caller_config[MitochondrialCallingKeys.INPUT_BAM])
        if not input_bam.endswith(".bam"):
            raise ValueError("GATK mitochondrial calling currently accepts BAM input only.")
        input_index = caller_config.get(MitochondrialCallingKeys.INPUT_BAM_INDEX)
        if input_index is not None:
            input_index = os.path.abspath(input_index)
            if not os.path.isfile(input_index):
                raise FileNotFoundError(f"BAM index does not exist: {input_index}")
        if not os.path.isfile(input_bam):
            raise FileNotFoundError(f"BAM does not exist: {input_bam}")

        def output(key: str) -> str:
            return cls._path(caller_config[key], output_dir)

        subset_bam = output(MitochondrialCallingKeys.SUBSET_BAM)
        unmapped_bam = output(MitochondrialCallingKeys.UNMAPPED_BAM)
        standard_bam = output(MitochondrialCallingKeys.STANDARD_ALIGNED_BAM)
        shifted_bam = output(MitochondrialCallingKeys.SHIFTED_ALIGNED_BAM)
        standard_raw = output(MitochondrialCallingKeys.STANDARD_RAW_VCF)
        shifted_raw = output(MitochondrialCallingKeys.SHIFTED_RAW_VCF)
        merged_raw = output(MitochondrialCallingKeys.MERGED_RAW_VCF)
        merged_stats = output(MitochondrialCallingKeys.MERGED_STATS)
        filtered_vcf = output(MitochondrialCallingKeys.FILTERED_VCF)
        split_vcf = output(MitochondrialCallingKeys.SPLIT_VCF)
        coverage_metrics = output(MitochondrialCallingKeys.COVERAGE_METRICS)
        duplicate_metrics = output(MitochondrialCallingKeys.DUPLICATE_METRICS)
        contamination_report = output(MitochondrialCallingKeys.CONTAMINATION_REPORT)
        contig = caller_config[MitochondrialCallingKeys.CONTIG_NAME]
        sample_name = Path(filtered_vcf).name.split(".")[0]
        read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"
        bind_paths = [{"source": str(Path(input_bam).parent), "target": str(Path(input_bam).parent), "mode": "ro"}]
        if input_index is not None:
            bind_paths.append(
                {"source": str(Path(input_index).parent), "target": str(Path(input_index).parent), "mode": "ro"}
            )
        runner = DockerRunner(device="cpu")

        cls._run(
            runner,
            cls._quote([
                "gatk", "PrintReads", "-I", input_bam, "-L", contig,
                "--read-filter", "MateOnSameContigOrNoMappedMateReadFilter",
                "--read-filter", "MateUnmappedAndUnmappedReadFilter",
            ] + (["--read-index", input_index] if input_index else []) + ["-O", subset_bam]), output_dir, bind_paths,
        )
        cls._run(
            runner,
            cls._quote([
                "gatk", "RevertSam", "-I", subset_bam, "-O", unmapped_bam,
                "--OUTPUT_BY_READGROUP", "false", "--SORT_ORDER", "queryname",
                "--RESTORE_ORIGINAL_QUALITIES", "false", "--ATTRIBUTE_TO_CLEAR", "FT",
                "--ATTRIBUTE_TO_CLEAR", "CO", "--VALIDATION_STRINGENCY", "LENIENT",
            ]), output_dir, bind_paths,
        )

        for reference, aligned_bam, metrics in (
            (refs["standard FASTA"], standard_bam, duplicate_metrics),
            (refs["shifted FASTA"], shifted_bam, f"{duplicate_metrics}.shifted"),
        ):
            temporary_bam = f"{aligned_bam}.bwa.bam"
            command = (
                cls._quote(["samtools", "fastq", "-n", "-0", "/dev/null", "-s", "/dev/null", unmapped_bam])
                + f" | bwa mem -p -R {shlex.quote(read_group)} {shlex.quote(reference)} /dev/stdin | "
                f"samtools sort -o {shlex.quote(temporary_bam)} -; "
                + cls._quote([
                    "gatk", "MarkDuplicates", "--INPUT", temporary_bam, "--OUTPUT", aligned_bam,
                    "--METRICS_FILE", metrics, "--CREATE_INDEX", "true",
                    "--VALIDATION_STRINGENCY", "SILENT",
                ])
            )
            cls._run(runner, command, output_dir, bind_paths)

        cls._run(
            runner,
            cls._quote([
                "gatk", "CollectWgsMetrics", "-I", standard_bam, "-R", refs["standard FASTA"],
                "-O", coverage_metrics, "--COVERAGE_CAP", "100000", "--VALIDATION_STRINGENCY", "SILENT",
            ]), output_dir, bind_paths,
        )
        coverage_option = []
        if caller_config.get(MitochondrialCallingKeys.MEDIAN_AUTOSOMAL_COVERAGE) is not None:
            coverage_option = ["--median-autosomal-coverage", caller_config[MitochondrialCallingKeys.MEDIAN_AUTOSOMAL_COVERAGE]]
        for reference, bam, interval, raw_vcf in (
            (refs["standard FASTA"], standard_bam, f"{contig}:576-16024", standard_raw),
            (refs["shifted FASTA"], shifted_bam, f"{contig}:8025-9144", shifted_raw),
        ):
            cls._run(
                runner,
                cls._quote([
                    "gatk", "Mutect2", "-R", reference, "-I", bam, "-L", interval,
                    "-O", raw_vcf, "--mitochondria-mode", "--max-reads-per-alignment-start", "75",
                    "--max-mnp-distance", "0", "--annotation", "StrandBiasBySample",
                ] + coverage_option), output_dir, bind_paths,
            )

        lifted_vcf = f"{shifted_raw}.lifted.vcf.gz"
        cls._run(
            runner,
            cls._quote([
                "gatk", "LiftoverVcf", "-I", shifted_raw, "-O", lifted_vcf,
                "-R", refs["standard FASTA"], "--CHAIN", refs["shift-back chain"],
                "--REJECT", f"{lifted_vcf}.rejected.vcf.gz",
            ]) + "; " + cls._quote([
                "gatk", "MergeVcfs", "-I", standard_raw, "-I", lifted_vcf, "-O", merged_raw,
            ]) + "; " + cls._quote([
                "gatk", "MergeMutectStats", "--stats", f"{standard_raw}.stats",
                "--stats", f"{shifted_raw}.stats", "-O", merged_stats,
            ]), output_dir, bind_paths,
        )

        preliminary_vcf = f"{filtered_vcf}.preliminary.vcf.gz"
        mutect_filtered_vcf = f"{filtered_vcf}.mutect.vcf.gz"
        filter_options = ["--mitochondria-mode", "--stats", merged_stats]
        if caller_config.get(MitochondrialCallingKeys.MIN_ALLELE_FRACTION) is not None:
            filter_options.extend(["--min-allele-fraction", caller_config[MitochondrialCallingKeys.MIN_ALLELE_FRACTION]])
        if caller_config.get(MitochondrialCallingKeys.F_SCORE_BETA) is not None:
            filter_options.extend(["--f-score-beta", caller_config[MitochondrialCallingKeys.F_SCORE_BETA]])
        contamination = caller_config.get(MitochondrialCallingKeys.CONTAMINATION_ESTIMATE, 0)
        filter_options.extend(["--contamination-estimate", contamination])
        cls._run(
            runner,
            cls._quote(["gatk", "FilterMutectCalls", "-R", refs["standard FASTA"], "-V", merged_raw, "-O", preliminary_vcf] + filter_options),
            output_dir, bind_paths,
        )
        if caller_config[MitochondrialCallingKeys.RUN_HAPLOCHECKER]:
            haplochecker_vcf = f"{filtered_vcf}.haplochecker.vcf.gz"
            cls._run(
                runner,
                cls._quote([
                    "gatk", "SelectVariants", "-V", preliminary_vcf, "--exclude-filtered",
                    "-O", haplochecker_vcf,
                ]), output_dir, bind_paths,
            )
            cls._run_haplochecker(output_dir, haplochecker_vcf, contamination_report, bind_paths)

        final_filter = cls._quote(
            [
                "gatk",
                "FilterMutectCalls",
                "-R",
                refs["standard FASTA"],
                "-V",
                merged_raw,
                "-O",
                mutect_filtered_vcf,
            ]
            + filter_options
        )
        if (
            caller_config[MitochondrialCallingKeys.RUN_HAPLOCHECKER]
            and MitochondrialCallingKeys.CONTAMINATION_ESTIMATE not in caller_config
        ):
            # Haplochecker reports heteroplasmy per major/minor haplogroup. This
            # matches the contamination calculation used by GATK's reference WDL.
            report_contamination = (
                "contamination=$(awk -F '\\t' 'NR == 1 { "
                "for (i = 1; i <= NF; i++) { if ($i == \"Contamination\") c=i; "
                "if ($i == \"MeanHetLevelMajor\") major=i; "
                "if ($i == \"MeanHetLevelMinor\") minor=i } next } "
                "NR == 2 { if ($c == \"YES\") { if ($(major) == 0) print $(minor); "
                "else print 1 - $(major) } else print 0; exit }' "
                f"{shlex.quote(contamination_report)}); "
                "[[ $contamination =~ ^(0|1)(\\.[0-9]+)?$ ]] || contamination=0; "
                + cls._quote(
                    [
                        "gatk",
                        "FilterMutectCalls",
                        "-R",
                        refs["standard FASTA"],
                        "-V",
                        merged_raw,
                        "-O",
                        mutect_filtered_vcf,
                    ]
                    + filter_options[:-2]
                )
                + ' --contamination-estimate "$contamination"'
            )
            final_filter = report_contamination

        cls._run(
            runner,
            final_filter + "; " + cls._quote([
                "gatk", "VariantFiltration", "-R", refs["standard FASTA"], "-V", mutect_filtered_vcf,
                "-O", filtered_vcf, "--apply-allele-specific-filters", "--mask", refs["blacklist"],
                "--mask-name", "blacklisted_site",
            ]) + "; " + cls._quote([
                "gatk", "LeftAlignAndTrimVariants", "-R", refs["standard FASTA"], "-V", filtered_vcf,
                "-O", split_vcf, "--split-multi-allelics", "--dont-trim-alleles", "--keep-original-ac",
            ]), output_dir, bind_paths,
        )
