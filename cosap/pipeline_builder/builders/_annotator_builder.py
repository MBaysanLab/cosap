from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import AnnotatorKeys, PipelineKeys
from ..._utils import join_paths
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from ._variantcaller_builder import VariantCaller


@dataclass
class Annotator(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    library: str
    sample_name: str = None
    name: str = None
    key: str = PipelineKeys.ANNOTATION
    next_step: _PipelineStep = None

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.input_step.name}_{self.library}"

        #Retrieve sample name from input step if possible.
        if self.sample_name is None:
            if self.input_step.__class__ == VariantCaller:
                self.sample_name = (
                    self.input_step.tumor
                    if self.input_step.tumor is not None
                    else self.input_step.normal
                )
            elif self.input_step.__class__ == Annotator:
                self.sample_name = self.input_step.sample_name
            else:
                raise Exception(
                    "Sample name cannot be read from input, please specify it by setting"
                    "sample_name argument."
                )

        self.input_step.next_step = self

    def _create_output_filename(self) -> str:
        if self.library.lower() in ["annovar", "intervar","cancervar"]:
            return FileFormats.ANNOTATION_OUTPUT.format(
                identification=self.name, custom_ext="txt"
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
                AnnotatorKeys.OUTPUT: join_paths(
                    OutputFolders.ANNOTATION, self.library, output_filename
                ),
                AnnotatorKeys.OUTPUT_DIR: OutputFolders.ANNOTATION,
            }
        }
        if self.library.lower() == "annovar":
            config[self.name][AnnotatorKeys.AVOUTPUT] = av_output_filename

        return config

    def get_output(self) -> str:
        config = self.get_config()
        return config[PipelineKeys.ANNOTATION][self.name][AnnotatorKeys.OUTPUT]

    def get_config(self) -> Dict:
        annotation_config = self._create_config()
        config = {PipelineKeys.ANNOTATION: annotation_config}
        return config
