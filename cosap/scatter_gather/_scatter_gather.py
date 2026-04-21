import glob
import os
import shutil
import subprocess
from contextlib import contextmanager
from itertools import chain, repeat
from subprocess import run

import shortuuid

from .._config import AppConfig
from .._pipeline_config import PipelineBaseKeys, PipelineKeys, VariantCallingKeys
from ..pipeline_builder import VariantCaller
from .utils import get_region_file_list, split_bam_by_intervals


class ScatterGather:
    @staticmethod
    def split_variantcaller_configs(
        config: dict, bed_file=None, split_bams: bool = False
    ) -> list[dict]:
        # If the number of threads is not suitable for parellelization, return the original config
        app_config = AppConfig()
        threads = app_config.MAX_THREADS_PER_JOB

        if int(threads) == 1:
            return [config]

        germline_bam = (
            config[VariantCallingKeys.GERMLINE_INPUT]
            if VariantCallingKeys.GERMLINE_INPUT in config.keys()
            else None
        )
        tumor_bam = (
            config[VariantCallingKeys.TUMOR_INPUT]
            if VariantCallingKeys.TUMOR_INPUT in config.keys()
            else None
        )

        if split_bams:
            if not germline_bam and not tumor_bam:
                raise ValueError(
                    "split_bams=True but neither germline nor tumor BAM is provided"
                )
            splitted_germline_bams = (
                split_bam_by_intervals(germline_bam, bed_file=bed_file)
                if germline_bam
                else repeat(None)
            )
            splitted_tumor_bams = (
                split_bam_by_intervals(tumor_bam, bed_file=bed_file)
                if tumor_bam
                else repeat(None)
            )
            bam_pairs = list(zip(splitted_germline_bams, splitted_tumor_bams))

        interval_files = get_region_file_list(
            file_type="interval_list",
            bed_file=bed_file,
            workdir=config[PipelineBaseKeys.OUTPUT_DIR],
        )
        splitted_configs = []

        for i in range(len(interval_files)):
            tmp_name = f"tmp{shortuuid.uuid()}"
            variant_caller = VariantCaller(
                library=config[VariantCallingKeys.LIBRARY],
                name=tmp_name,
                bed_file=interval_files[i],
                params=config[PipelineBaseKeys.PARAMS],
                gvcf=config[VariantCallingKeys.OUTPUT_TYPE] == "GVCF",
                output_dir=config[PipelineBaseKeys.OUTPUT_DIR],
            )
            cfg = variant_caller.get_config()[PipelineKeys.VARIANT_CALLING][tmp_name]

            if VariantCallingKeys.GERMLINE_INPUT in config.keys():
                cfg[VariantCallingKeys.GERMLINE_INPUT] = (
                    bam_pairs[i][0] if split_bams else germline_bam
                )
            if VariantCallingKeys.TUMOR_INPUT in config.keys():
                cfg[VariantCallingKeys.TUMOR_INPUT] = (
                    bam_pairs[i][1] if split_bams else tumor_bam
                )
            cfg = dict(config, **cfg)
            splitted_configs.append(cfg)

        return splitted_configs

    @staticmethod
    @contextmanager
    def cleanup_on_exit(filepath):
        try:
            yield
        except KeyboardInterrupt:
            raise
        except Exception:
            if os.path.exists(filepath):
                os.remove(filepath)
            # Also clean up associated index files
            for ext in [".tbi", ".idx"]:
                idx_path = filepath + ext
                if os.path.exists(idx_path):
                    os.remove(idx_path)
            raise

    @staticmethod
    def gather_vcfs(configs: list, output_path, mode="vcf", cwd=None):
        GVCF_MODE = "gvcf"
        VCF_MODE = "vcf"

        if mode.lower() == VCF_MODE:
            vcfs = [
                cfg[VariantCallingKeys.UNFILTERED_VARIANTS_OUTPUT] for cfg in configs
            ]
            command = ["gatk", "MergeVcfs", "-O", output_path]
        elif mode.lower() == GVCF_MODE:
            vcfs = [cfg[VariantCallingKeys.GVCF_OUTPUT] for cfg in configs]
            command = ["gatk", "SortVcf", "-O", output_path]
        else:
            raise ValueError(f"Mode {mode} not supported")

        if not vcfs:
            raise RuntimeError(
                "No VCF files to gather. The scatter step may have produced no output. "
                "Check that interval files exist and scatter jobs completed successfully."
            )

        # Verify all input VCF files exist
        missing_vcfs = [
            vcf for vcf in vcfs if not os.path.isfile(os.path.join(cwd, vcf) if cwd else vcf)
        ]
        if missing_vcfs:
            raise RuntimeError(
                f"Missing VCF files for gather step: {missing_vcfs}. "
                "Scatter jobs may have failed to produce output."
            )

        command.extend(list(chain(*zip(repeat("-I"), vcfs))))

        with ScatterGather.cleanup_on_exit(output_path):
            try:
                run(command, cwd=cwd, check=True)
            except subprocess.CalledProcessError as e:
                print(e)
                raise RuntimeError("Failed to gather VCFs") from e

        # If we reach here, the command was successful, output file is preserved

    @staticmethod
    def clean_temp_files(path):
        # Scatter shards use identifications like "tmp<shortuuid>", which can
        # appear either at the start of a filename (e.g. "tmpABC.g.vcf.gz")
        # or embedded after a fixed prefix (e.g. "all_unfiltered_tmpABC.vcf.gz").
        # Match both, plus their associated index files.
        temp_files = set(glob.glob(f"{path}/tmp*")) | set(glob.glob(f"{path}/*_tmp*"))
        for tmp in temp_files:
            if os.path.isfile(tmp):
                os.remove(tmp)
            elif os.path.isdir(tmp):
                shutil.rmtree(tmp)

    @staticmethod
    def split_bam_process_configs(config: dict) -> list[dict]:
        splitted_configs = []

        input_bam = config[PipelineBaseKeys.INPUT]
        splitted_bams = split_bam_by_intervals(input_bam)
        for bam in splitted_bams:
            cnf = {}
            cnf[PipelineBaseKeys.INPUT] = bam
            splitted_configs.append(cnf)

        return splitted_configs
