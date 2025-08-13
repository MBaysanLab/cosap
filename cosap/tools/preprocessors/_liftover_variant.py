import os
from subprocess import run
from typing import Dict

from ..._library_paths import LibraryPaths
from ..._pipeline_config import LiftoverKeys
from ..._utils import convert_to_absolute_path
from ._preprocessors import _PreProcessable, _Preprocessor


class LiftoverVariant(_Preprocessor, _PreProcessable):
    @classmethod
    def _create_command(cls, liftover_config: Dict, library_paths: LibraryPaths):
        """
        Create GATK LiftoverVariant command.
        """
        command = [
            "gatk",
            "LiftoverVcf",
            "-I",
            convert_to_absolute_path(liftover_config[LiftoverKeys.INPUT]),
            "-O",
            liftover_config[LiftoverKeys.OUTPUT],
            "--CHAIN",
            library_paths.CHAIN_FILE,
            "--REFERENCE_SEQUENCE",
            library_paths.REF_FASTA,
            "--REJECT",
            liftover_config[LiftoverKeys.REJECT_OUTPUT],
        ]

        return command

    @classmethod
    def _run_liftover(cls, liftover_config: Dict, library_paths: LibraryPaths):
        """
        Execute the GATK LiftoverVariant command.
        """
        # Create output directory if it doesn't exist
        output_dir = os.path.abspath(liftover_config[LiftoverKeys.OUTPUT_DIR])
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

        command = cls._create_command(
            liftover_config=liftover_config, library_paths=library_paths
        )
        run(command, cwd=liftover_config[LiftoverKeys.OUTPUT_DIR])

    @classmethod
    def run_preprocessor(cls, liftover_config: Dict, device="cpu"):
        """
        Main entry point to run the LiftoverVariant tool.

        Args:
            liftover_config: Configuration dictionary containing input/output paths
            device: Execution device - "cpu" or "gpu" (currently only CPU supported)
        """
        library_paths = LibraryPaths()

        if device == "cpu":
            cls._run_liftover(
                liftover_config=liftover_config, library_paths=library_paths
            )
        else:
            raise NotImplementedError(
                f"GPU execution not yet supported for LiftoverVariant. "
                f"Requested device: {device}"
            )
