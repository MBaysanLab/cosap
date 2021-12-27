from typing import Dict

from .pipeline_runner import PipelineRunner


class Cosap:
    def run_pipeline(self, pipeline_config: Dict, snakemake: bool):
        PipelineRunner.run_pipeline(pipeline_config=pipeline_config, snakemake=snakemake)

    def mock_run_pipeline(self, pipeline_config: Dict):
        # TODO: should this be here?
        raise NotImplementedError()
