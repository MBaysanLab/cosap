from dataclasses import dataclass
from typing import Dict, List

from ..._formats import FileFormats, OutputFolders
from ..._pipeline_config import AnnotatorKeys, PipelineKeys
from ._pipeline_steps import _IPipelineStep, _PipelineStep
from ..._utils import join_paths


@dataclass
class Annotator(_IPipelineStep, _PipelineStep):
    input_step: _PipelineStep
    library: str
    name: str = None
    key: str = PipelineKeys.ANNOTATION

    def __post_init__(self):
        if self.name is None:
            self.name = f"{self.input_step.name}_{self.library}"

    def _create_config(self) -> Dict:
        output_filename = FileFormats.ANNOTATING_OUTPUT.format(identification=self.name)
        av_output_filename = FileFormats.ANNOVAR_OUTPUT.format(identification=self.name)
        config = {
            self.name: {
                AnnotatorKeys.LIBRARY: self.library,
                AnnotatorKeys.INPUT: self.input_step.get_output(),
                AnnotatorKeys.OUTPUT: output_filename,
                AnnotatorKeys.OUTPUT_DIR: OutputFolders.ANNOTATION,
            }
        }
        if self.library.lower() == "annovar":
            config[self.name][AnnotatorKeys.AVOUTPUT] = av_output_filename

        return config

    def get_output(self) -> str:
        config = self.get_config()
        return join_paths(
            config[PipelineKeys.ANNOTATION][self.name][AnnotatorKeys.OUTPUT_DIR],
            config[PipelineKeys.ANNOTATION][self.name][AnnotatorKeys.OUTPUT],
            )

    def get_config(self) -> Dict:
        annotation_config = self._create_config()
        config = {PipelineKeys.ANNOTATION: annotation_config}
        return config
