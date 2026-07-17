import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from ..._docker_images import DockerImages
from ..._pipeline_config import CNVCallingKeys
from ...runners.docker_runner import DockerRunner
from ._cnv_callers import _CNVCaller


class FREEC(_CNVCaller):
    """Control-FREEC WES CNA caller.

    Supports matched sample/control WES and control-free WES. BAF/LOH output is
    intentionally not wired in this first implementation.
    """

    DEFAULT_READ_COUNT_THRESHOLD = 50

    @staticmethod
    def _as_abs_path(path: str, field_name: str) -> str:
        if not path:
            raise ValueError(f"FREEC requires '{field_name}'.")
        return str(Path(path).expanduser().resolve())

    @staticmethod
    def _optional_abs_path(path: Optional[str]) -> Optional[str]:
        if not path:
            return None
        return str(Path(path).expanduser().resolve())

    @staticmethod
    def _sample_stem(path: str) -> str:
        name = Path(path).name
        for suffix in (".bam", ".sam", ".pileup.gz", ".pileup", ".gz"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return Path(name).stem

    @classmethod
    def _normalise_input_list(cls, value) -> List[str]:
        if value is None:
            return []
        if isinstance(value, (str, os.PathLike)):
            return [str(value)]
        return [str(item) for item in value]

    @classmethod
    def _validate_inputs(cls, caller_config: Dict) -> Dict:
        params = dict(caller_config.get(CNVCallingKeys.PARAMS, {}))
        tumor_inputs = cls._normalise_input_list(
            caller_config.get(CNVCallingKeys.TUMOR_INPUT)
        )
        normal_inputs = cls._normalise_input_list(
            caller_config.get(CNVCallingKeys.NORMAL_INPUT)
        )

        if len(tumor_inputs) != 1:
            raise ValueError("FREEC currently requires exactly one sample/tumor input.")
        if len(normal_inputs) > 1:
            raise ValueError(
                "FREEC currently supports at most one matched control input."
            )

        sample_input = cls._as_abs_path(tumor_inputs[0], CNVCallingKeys.TUMOR_INPUT)
        control_input = (
            cls._optional_abs_path(normal_inputs[0]) if normal_inputs else None
        )
        capture_regions = cls._as_abs_path(
            params.get("capture_regions") or caller_config.get(CNVCallingKeys.BED_FILE),
            "capture_regions",
        )
        chr_len_file = cls._as_abs_path(params.get("chr_len_file"), "chr_len_file")

        if sample_input.endswith(".cram") or (
            control_input and control_input.endswith(".cram")
        ):
            raise ValueError(
                "FREEC does not document CRAM input support; provide BAM, SAM, or pileup."
            )

        gc_content_profile = cls._optional_abs_path(params.get("gc_content_profile"))
        chr_files = cls._optional_abs_path(params.get("chr_files"))
        if not control_input and not (gc_content_profile or chr_files):
            raise ValueError(
                "FREEC control-free WES requires either params.gc_content_profile or params.chr_files."
            )

        return {
            "params": params,
            "sample_input": sample_input,
            "control_input": control_input,
            "capture_regions": capture_regions,
            "chr_len_file": chr_len_file,
            "gc_content_profile": gc_content_profile,
            "chr_files": chr_files,
        }

    @classmethod
    def _config_lines(cls, resolved: Dict, output_dir: str) -> List[str]:
        params = resolved["params"]
        lines = [
            "[general]",
            f"chrLenFile = {resolved['chr_len_file']}",
            "window = 0",
            f"ploidy = {params.get('ploidy', 2)}",
            f"outputDir = {output_dir}",
            f"maxThreads = {int(params.get('threads', params.get('max_threads', 1)))}",
            (
                "noisyData = TRUE"
                if params.get("noisy_data", True)
                else "noisyData = FALSE"
            ),
            (
                "printNA = FALSE"
                if not params.get("print_na", False)
                else "printNA = TRUE"
            ),
            f"readCountThreshold = {int(params.get('read_count_threshold', cls.DEFAULT_READ_COUNT_THRESHOLD))}",
        ]

        optional_general = {
            "break_point_type": "breakPointType",
            "break_point_threshold": "breakPointThreshold",
            "minimal_subclone_presence": "minimalSubclonePresence",
            "sex": "sex",
            "degree": "degree",
            "force_gc_content_normalization": "forceGCcontentNormalization",
        }
        for key, freec_key in optional_general.items():
            if key in params and params[key] is not None:
                lines.append(f"{freec_key} = {params[key]}")

        if resolved["gc_content_profile"]:
            lines.append(f"GCcontentProfile = {resolved['gc_content_profile']}")
        if resolved["chr_files"]:
            lines.append(f"chrFiles = {resolved['chr_files']}")
        if params.get("gem_mappability_file"):
            lines.append(
                f"gemMappabilityFile = {cls._as_abs_path(params['gem_mappability_file'], 'gem_mappability_file')}"
            )
            lines.append("uniqueMatch = TRUE")

        input_format = params.get("input_format", "BAM")
        mate_orientation = params.get("mate_orientation", 0)
        lines.extend(
            [
                "",
                "[sample]",
                f"mateFile = {resolved['sample_input']}",
                f"inputFormat = {input_format}",
                f"mateOrientation = {mate_orientation}",
            ]
        )

        if resolved["control_input"]:
            lines.extend(
                [
                    "",
                    "[control]",
                    f"mateFile = {resolved['control_input']}",
                    f"inputFormat = {params.get('control_input_format', input_format)}",
                    f"mateOrientation = {params.get('control_mate_orientation', mate_orientation)}",
                ]
            )

        lines.extend(
            [
                "",
                "[target]",
                f"captureRegions = {resolved['capture_regions']}",
                "",
            ]
        )
        return lines

    @classmethod
    def _write_config(cls, caller_config: Dict, resolved: Dict, output_dir: str) -> str:
        config_path = caller_config.get(CNVCallingKeys.FREEC_CONFIG) or os.path.join(
            output_dir, "freec_config.txt"
        )
        with open(config_path, "w") as handle:
            handle.write("\n".join(cls._config_lines(resolved, output_dir)))
        caller_config[CNVCallingKeys.FREEC_CONFIG] = config_path
        return config_path

    @classmethod
    def _bind_paths(
        cls, paths: Iterable[Optional[str]], output_dir: str
    ) -> List[Dict[str, str]]:
        bind_dirs = {output_dir}
        for path in paths:
            if not path:
                continue
            resolved = Path(path).resolve()
            bind_dirs.add(str(resolved if resolved.is_dir() else resolved.parent))
            if resolved.suffix == ".bam":
                for index_suffix in (".bai", ".csi"):
                    index = Path(str(resolved) + index_suffix)
                    if index.exists():
                        bind_dirs.add(str(index.parent.resolve()))
        return [
            {
                "source": bind_dir,
                "target": bind_dir,
                "mode": "rw" if bind_dir == output_dir else "ro",
            }
            for bind_dir in sorted(bind_dirs)
        ]

    @staticmethod
    def _parse_info(path: str) -> Dict[str, str]:
        values = {}
        if not path or not os.path.exists(path):
            return values
        with open(path) as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                if "=" in stripped:
                    key, value = stripped.split("=", 1)
                elif "\t" in stripped:
                    key, value = stripped.split("\t", 1)
                else:
                    continue
                values[key.strip()] = value.strip()
        return values

    @classmethod
    def _is_control_free_run(cls, caller_config: Dict) -> bool:
        return not cls._normalise_input_list(
            caller_config.get(CNVCallingKeys.NORMAL_INPUT)
        )

    @classmethod
    def _find_output_file(cls, output_dir: str, suffix: str) -> Optional[str]:
        candidates = sorted(Path(output_dir).glob(f"*{suffix}"))
        return str(candidates[0]) if candidates else None

    @classmethod
    def _validate_run_outputs(cls, caller_config: Dict, output_dir: str) -> str:
        info_path = cls._find_output_file(output_dir, "_info.txt")
        raw_cnvs = cls._find_output_file(output_dir, "_CNVs")
        if not raw_cnvs:
            raise RuntimeError(
                f"FREEC completed but no *_CNVs file was found in {output_dir}."
            )
        if info_path:
            info = cls._parse_info(info_path)
            if info.get("Good_Polynomial_Fit", "True").lower() == "false":
                raise RuntimeError(
                    f"FREEC reported Good_Polynomial_Fit=False in {info_path}."
                )
            if (
                cls._is_control_free_run(caller_config)
                and info.get("CGcontent_Used", "True").lower() != "true"
            ):
                raise RuntimeError(
                    "FREEC reported CGcontent_Used=False for a control-free WES run. "
                    "Control-FREEC v11.6 disables normalization for targeted data "
                    "without a control sample, so this output is not accepted. "
                    "Provide a matched control input or use a caller that supports "
                    "control-free targeted normalization."
                )
            caller_config[CNVCallingKeys.FREEC_INFO] = info_path
        caller_config[CNVCallingKeys.FREEC_RAW_CNVS] = raw_cnvs
        return raw_cnvs

    @classmethod
    def _normalise_cnvs(cls, raw_cnvs: str, output_path: str) -> None:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(raw_cnvs) as in_handle, open(output_path, "w") as out_handle:
            out_handle.write("chrom\tstart\tend\tcopy_number\tstatus\textra\n")
            for line in in_handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split("\t")
                if len(parts) < 5:
                    parts = stripped.split()
                if len(parts) < 5:
                    continue
                chrom, start, end, copy_number, status = parts[:5]
                extra = "|".join(parts[5:])
                out_handle.write(
                    f"{chrom}\t{start}\t{end}\t{copy_number}\t{status}\t{extra}\n"
                )

    @classmethod
    def call(cls, caller_config: dict, device: str = "cpu", max_workers=None):
        del max_workers
        output_dir = cls._as_abs_path(
            caller_config[CNVCallingKeys.OUTPUT_DIR], CNVCallingKeys.OUTPUT_DIR
        )
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, caller_config[CNVCallingKeys.OUTPUT])

        resolved = cls._validate_inputs(caller_config)
        config_path = cls._write_config(caller_config, resolved, output_dir)
        bind_paths = cls._bind_paths(
            [
                config_path,
                resolved["sample_input"],
                resolved["control_input"],
                resolved["capture_regions"],
                resolved["chr_len_file"],
                resolved["gc_content_profile"],
                resolved["chr_files"],
                resolved["params"].get("gem_mappability_file"),
            ],
            output_dir,
        )

        runner = DockerRunner(device=device)
        return_code = runner.run(
            DockerImages.FREEC,
            ["freec", "-conf", config_path],
            workdir=output_dir,
            paths_to_bind=bind_paths,
        )
        if return_code not in (0, None):
            raise RuntimeError(f"FREEC failed with return code {return_code}.")

        raw_cnvs = cls._validate_run_outputs(caller_config, output_dir)
        cls._normalise_cnvs(raw_cnvs, output_path)
        return output_path
