from .pipeline_runner import PipelineRunner


class Cosap:
    def run_pipeline(self, pipeline_config: Dict):
        PipelineRunner.run_pipeline(pipeline_config=pipeline_config)

    def mock_run_pipeline(self, pipeline_config: Dict):
        # TODO: should this be here?
        pass
