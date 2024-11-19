from typing import Dict


from .runners import PipelineRunner

class Cosap:
    def run_pipeline(self, pipeline_config: Dict, snakemake: str):
        PipelineRunner.run_pipeline(pipeline_config=pipeline_config, backend=snakemake)

    def mock_run_pipeline(self, pipeline_config: Dict):
        # TODO: should this be here?
        raise NotImplementedError()

