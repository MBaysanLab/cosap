from typing import Dict

from .pipeline_runner import PipelineRunner


class Cosap:
    def run_pipeline(self, pipeline_config: Dict):
        PipelineRunner.run_pipeline(pipeline_config=pipeline_config)

    def run_pipeline_snakemake(self, pipeline_config: Dict, workdir: str):
        PipelineRunner.run_pipeline_snakemake(
            pipeline_config=pipeline_config, workdir=workdir
        )

    def mock_run_pipeline(self, pipeline_config: Dict):
        # TODO: should this be here?
        raise NotImplementedError()
