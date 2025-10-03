from dataclasses import dataclass
from typing import Dict, Optional

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import AnnotatorKeys, PipelineKeys, VariantCallingKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from ._variantcaller_builder import VariantCaller


@dataclass
class Annotator(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    library: str
    sample_name: Optional[str] = None
    name: Optional[str] = None
    key: str = PipelineKeys.ANNOTATION
    next_step: Optional[_PipelineStep] = None
    input_type: str = "vcf"
    params: Optional[Dict] = None
    output_dir: Optional[str] = None

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.input_step.name}_{self.library}"

        self.library = self.library.upper()

        # Retrieve sample name from input step if possible.
        if self.sample_name is None:
            if self.input_step.__class__ == VariantCaller:
                self.sample_name = (
                    self.input_step.params[VariantCallingKeys.TUMOR_SAMPLE_NAME]
                    if self.input_step.tumor is not None
                    else self.input_step.params[VariantCallingKeys.GERMLINE_SAMPLE_NAME]
                )
            elif self.input_step.__class__ == Annotator:
                self.sample_name = self.input_step.sample_name
            else:
                raise Exception(
                    "Sample name cannot be read from input, please specify it by setting"
                    "sample_name argument."
                )

        self.input_step.next_step = self

        if self.output_dir is None:
            self.output_dir = join_paths(OutputFolders.ANNOTATION, self.library)

    def _create_output_filename(self) -> str:
        if self.library.lower() in ["annovar", "intervar", "cancervar"]:
            return FileFormats.ANNOTATION_OUTPUT.format(
                identification=self.name, custom_ext="txt"
            )
        elif self.library.lower() in ["annotsv"]:
            return FileFormats.ANNOTATION_OUTPUT.format(
                identification=self.name, custom_ext="tsv"
            )
        elif self.library.lower() in ["pharmcat", "exomiser"]:
            return FileFormats.ANNOTATION_OUTPUT.format(
                identification=self.name, custom_ext="json"
            )
        else:
            return FileFormats.ANNOTATION_OUTPUT.format(
                identification=self.name, custom_ext="vcf"
            )

    def _create_config(self) -> Dict:
        output_filename = self._create_output_filename()
        av_output_filename = FileFormats.ANNOVAR_OUTPUT.format(
            identification=self.name, sample=self.sample_name
        )
        config = {
            self.name: {
                AnnotatorKeys.LIBRARY: self.library,
                AnnotatorKeys.INPUT: self.input_step.get_output(),
                AnnotatorKeys.OUTPUT: output_filename,
                AnnotatorKeys.OUTPUT_DIR: self.output_dir,
                AnnotatorKeys.INPUT_TYPE: self.input_type,
                AnnotatorKeys.LOG_FILE: self.log_file,
                AnnotatorKeys.PARAMS: self.params if self.params is not None else {},
            }
        }
        if self.library.lower() == "annovar":
            config[self.name][AnnotatorKeys.AVOUTPUT] = av_output_filename

        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(
            self.output_dir, config[self.key][self.name][AnnotatorKeys.OUTPUT]
        )

    def get_config(self, core_only: bool = False) -> Dict:
        annotation_config = self._create_config()
        if core_only:
            return annotation_config[self.name]
        config = {self.key: annotation_config}
        return config
