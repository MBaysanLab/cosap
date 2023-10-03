import os
from typing import Dict, List

import yaml

from .._config import AppConfig
from .._pipeline_config import MappingKeys, PipelineKeys, VariantCallingKeys
from .._utils import join_paths
from ..mappers import MapperFactory
from ..preprocessors import (BamIndexer, BamMerger, BaseRecalibrator,
                             MarkDuplicate, SamtoolsSorter)
from ..variant_callers import VariantCallerFactory
from ._snakemake_runner import SnakemakeRunner


class PipelineRunner:
    def validate_pipeline_config(self, pipeline_config: Dict):
        raise NotImplementedError()

    def map(self, mapping_configs: List):
        for config in mapping_configs:
            mapper = MapperFactory.create(config[MappingKeys.LIBRARY])
            mapper.map(config)

    def sort(self, sorting_configs: List):
        for config in sorting_configs:
            SamtoolsSorter.sort(config)

    def create_index(self, indexing_configs: List):
        # TODO: should check if each element is ready to be processed
        # maybe for all steps it can be done and even a mock run can be made this way
        for config in indexing_configs:
            BamIndexer.create_index(config)

    def merge(self, merge_config: List):
        # TODO: if there is a single bam skip this
        for config in merge_config:
            BamMerger.merge(config)

    def calibrate(self, calibrate_config: List):
        for config in calibrate_config:
            BaseRecalibrator.calibrate(config)

    def mark_duplicates(self, duplicates_config: List):
        for config in duplicates_config:
            MarkDuplicate.mark(config)

    def call_variants(self, variant_configs: List):
        for config in variant_configs:
            caller = VariantCallerFactory.create(config[VariantCallingKeys.LIBRARY])
            caller.call_variants(config)

    def _write_config_to_yaml(self, config_yaml_path, pipeline_config):
        with open(config_yaml_path, "w") as config_yaml:
            yaml.dump(pipeline_config, config_yaml, default_flow_style=False)

    def run_pipeline(self, pipeline_config: Dict, runner: str = "snakemake") -> str:
        workdir = pipeline_config[PipelineKeys.WORKDIR]
        config_yaml_path = os.path.normpath(
            join_paths(
                workdir, f"{pipeline_config[PipelineKeys.CREATION_DATE]}_config.yaml"
            )
        )
        self._write_config_to_yaml(config_yaml_path, pipeline_config)

        if runner.lower() == "snakemake":
            snakemake_runner = SnakemakeRunner(
                pipeline_config=config_yaml_path,
                workdir=pipeline_config[PipelineKeys.WORKDIR],
            )
            snakemake_runner.run_snakemake_pipeline()

        else:
            self.validate_pipeline_config(pipeline_config)

            # TODO: add trimming
            self.map(pipeline_config[PipelineKeys.MAPPING])
            self.sort(pipeline_config[PipelineKeys.SORTING])
            self.create_index(pipeline_config[PipelineKeys.INDEX])
            self.merge(pipeline_config[PipelineKeys.MERGE])
            self.calibrate(pipeline_config[PipelineKeys.CALIBRATE])

            self.call_variants(pipeline_config[PipelineKeys.VARIANT_CALLING])

        return pipeline_config
