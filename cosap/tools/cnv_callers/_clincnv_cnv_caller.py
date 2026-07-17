import csv
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from ..._docker_images import DockerImages
from ..._pipeline_config import CNVCallingKeys
from ...runners.docker_runner import DockerRunner
from ._cnv_callers import _CNVCaller


class ClinCNV(_CNVCaller):
    """ClinCNV germline caller operating on a cohort coverage matrix."""

    DEFAULT_CONTAINER_SCRIPT_DIR = "/opt/ClinCNV"

    @staticmethod
    def _as_abs_path(path: str, field_name: str) -> str:
        if not path:
            raise ValueError(f"ClinCNV requires '{field_name}'.")
        resolved = Path(path).expanduser().resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"ClinCNV {field_name} does not exist: {resolved}")
        return str(resolved)

    @classmethod
    def _validate_inputs(cls, caller_config: Dict) -> Dict:
        params = dict(caller_config.get(CNVCallingKeys.PARAMS, {}))
        coverage_matrix = cls._as_abs_path(
            params.get("coverage_matrix"), "params.coverage_matrix"
        )
        annotated_bed = cls._as_abs_path(
            params.get("annotated_bed") or caller_config.get(CNVCallingKeys.BED_FILE),
            "annotated_bed",
        )
        host_script_dir = None
        script_dir = params.get("script_dir") or cls.DEFAULT_CONTAINER_SCRIPT_DIR
        if params.get("script_dir"):
            host_script_dir = cls._as_abs_path(
                params.get("script_dir"), "params.script_dir"
            )
            script = Path(host_script_dir) / "clinCNV.R"
            if not script.exists():
                raise FileNotFoundError(f"ClinCNV script was not found: {script}")
            script_dir = host_script_dir
        sample_name = params.get("sample_name")
        if not sample_name:
            raise ValueError("ClinCNV requires 'params.sample_name'.")

        with open(coverage_matrix) as handle:
            header = handle.readline().rstrip("\n").split("\t")
        if sample_name not in header[3:]:
            raise ValueError(
                f"ClinCNV sample_name '{sample_name}' is not a coverage-matrix column."
            )

        return {
            "params": params,
            "coverage_matrix": coverage_matrix,
            "annotated_bed": annotated_bed,
            "script_dir": script_dir,
            "host_script_dir": host_script_dir,
            "sample_name": sample_name,
        }

    @classmethod
    def _command(cls, resolved: Dict, output_dir: str) -> List[str]:
        params = resolved["params"]
        command = [
            params.get("rscript", "Rscript"),
            str(Path(resolved["script_dir"]) / "clinCNV.R"),
            "--normal",
            resolved["coverage_matrix"],
            "--bed",
            resolved["annotated_bed"],
            "--out",
            output_dir,
            "--normalSample",
            resolved["sample_name"],
            "--folderWithScript",
            resolved["script_dir"],
        ]
        optional_params = {
            "score_g": "--scoreG",
            "length_g": "--lengthG",
            "max_num_germ_cnvs": "--maxNumGermCNVs",
            "max_num_iter": "--maxNumIter",
            "threads": "--numberOfThreads",
            "visualization_igv": "--visulizationIGV",
            "hg38": "--hg38",
        }
        for key, flag in optional_params.items():
            if key in params and params[key] is not None:
                command.extend([flag, str(params[key])])
        return command

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
        return [
            {
                "source": bind_dir,
                "target": bind_dir,
                "mode": "rw" if bind_dir == output_dir else "ro",
            }
            for bind_dir in sorted(bind_dirs)
        ]

    @staticmethod
    def _raw_output_path(output_dir: str, sample_name: str) -> str:
        return str(
            Path(output_dir) / "normal" / sample_name / f"{sample_name}_cnvs.tsv"
        )

    @classmethod
    def _normalise_cnvs(cls, raw_cnvs: str, output_path: str) -> None:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(raw_cnvs) as in_handle, open(output_path, "w") as out_handle:
            out_handle.write("chrom\tstart\tend\tcopy_number\tstatus\textra\n")
            rows = (line for line in in_handle if not line.startswith("##"))
            for row in csv.DictReader(rows, delimiter="\t"):
                copy_number = row["CN_change"]
                try:
                    copy_number_float = float(copy_number)
                except ValueError:
                    continue
                if copy_number_float == 2:
                    continue
                status = "loss" if copy_number_float < 2 else "gain"
                extra = "|".join(
                    f"{key}={value}"
                    for key, value in row.items()
                    if key not in {"#chr", "start", "end", "CN_change"}
                )
                out_handle.write(
                    f"{row['#chr']}\t{row['start']}\t{row['end']}\t"
                    f"{copy_number}\t{status}\t{extra}\n"
                )

    @classmethod
    def call(cls, caller_config: dict, device: str = "cpu", max_workers=None):
        del max_workers
        output_dir = str(
            Path(caller_config[CNVCallingKeys.OUTPUT_DIR]).expanduser().resolve()
        )
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, caller_config[CNVCallingKeys.OUTPUT])
        resolved = cls._validate_inputs(caller_config)
        bind_paths = cls._bind_paths(
            [
                resolved["coverage_matrix"],
                resolved["annotated_bed"],
                resolved["host_script_dir"],
            ],
            output_dir,
        )
        runner = DockerRunner(device=device)
        return_code = runner.run(
            resolved["params"].get("image") or DockerImages.CLINCNV,
            cls._command(resolved, output_dir),
            workdir=output_dir,
            paths_to_bind=bind_paths,
        )
        if return_code not in (0, None):
            raise RuntimeError(f"ClinCNV failed with return code {return_code}.")

        raw_cnvs = cls._raw_output_path(output_dir, resolved["sample_name"])
        if not os.path.exists(raw_cnvs):
            raise RuntimeError(
                f"ClinCNV completed but no call file was found: {raw_cnvs}"
            )
        caller_config[CNVCallingKeys.CLINCNV_RAW_CNVS] = raw_cnvs
        cls._normalise_cnvs(raw_cnvs, output_path)
        return output_path
