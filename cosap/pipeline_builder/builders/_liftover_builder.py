from dataclasses import dataclass, field
from typing import Dict, Optional

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import LiftoverKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep


@dataclass
class LiftoverVariant(_IPipelineStep, _PipelineStep):
    """
    Builder for GATK LiftoverVcf tool configuration.

    This builder creates configuration for lifting over variants from one reference
    genome to another using GATK's LiftoverVcf tool.

    Example usage:
        # Assuming you have a variant caller step
        variant_caller = VariantCaller(...)

        # Create liftover step to convert from hg19 to hg38
        liftover = LiftoverVariant(
            input_step=variant_caller,
            reject_output="rejected_variants"  # Optional
        )

        # Get the configuration
        config = liftover.get_config()

        # Get the output file path
        output_file = liftover.get_output()
    """

    input_step: _PipelineStep
    name: Optional[str] = None
    params: Dict = field(default_factory=dict)
    key: str = PipelineKeys.LIFTOVER
    next_step: Optional[_PipelineStep] = None
    output_dir: Optional[str] = None

    def __post_init__(self):
        self.key = PipelineKeys.LIFTOVER

        if self.name is None:
            self.name = f"{self.input_step.name}_liftover"

        self.input_step.next_step = self

        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.PREPROCESSOR, "liftover")

    def _create_config(self) -> Dict:
        """Create the configuration dictionary for liftover."""
        output_filename = FileFormats.LIFTOVER_OUTPUT.format(identification=self.name)
        reject_filename = FileFormats.LIFTOVER_REJECT_OUTPUT.format(
            identification=self.name
        )

        config = {
            self.name: {
                LiftoverKeys.INPUT: self.input_step.get_output(),
                LiftoverKeys.OUTPUT: output_filename,
                LiftoverKeys.OUTPUT_DIR: self.output_dir,
                LiftoverKeys.LOG_FILE: self.log_file,
                LiftoverKeys.PARAMS: self.params,
                LiftoverKeys.REJECT_OUTPUT: reject_filename,
            },
        }

        return config

    def get_output(self) -> str:
        """Get the output file path for the liftover step."""
        config = self.get_config()
        return join_paths(
            self.output_dir, config[self.key][self.name][LiftoverKeys.OUTPUT]
        )

    def get_reject_output(self) -> Optional[str]:
        """Get the reject output file path if specified."""
        if self.reject_output is not None:
            config = self.get_config()
            return join_paths(
                self.output_dir, config[self.key][self.name][LiftoverKeys.REJECT_OUTPUT]
            )
        return None

    def get_config(self, core_only: bool = False) -> Dict:
        """Get the complete configuration for this liftover step."""
        liftover_config = self._create_config()
        if core_only:
            return liftover_config[self.name]
        config = {self.key: liftover_config}
        return config
