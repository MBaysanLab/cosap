import os
from typing import Dict, List
import yaml
from subprocess import PIPE, Popen, check_output, run

from .._pipeline_config import MappingKeys, PipelineKeys, VariantCallingKeys
from ..mappers import MapperFactory
from ..preprocessors import (
    BamIndexer,
    BamMerger,
    MarkDuplicate,
    SamtoolsSorter,
    BaseRecalibrator,
)
from ..variant_callers import VariantCallerFactory
from .._config import AppConfig

from .._utils import join_paths


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

    def run_pipeline(self, pipeline_config: Dict):
        self.validate_pipeline_config(pipeline_config)

        # TODO: add trimming
        self.map(pipeline_config[PipelineKeys.MAPPING])
        self.sort(pipeline_config[PipelineKeys.SORTING])
        self.index(pipeline_config[PipelineKeys.INDEX])
        self.merge(pipeline_config[PipelineKeys.MERGE])
        self.calibrate(pipeline_config[PipelineKeys.CALIBRATE])

        self.call_variants(pipeline_config[PipelineKeys.VARIANT_CALLING])

    def run_pipeline_snakemake(self, pipeline_config: Dict, workdir: str):
        config = pipeline_config
        config[PipelineKeys.WORKDIR] = workdir

        config_yaml_path = join_paths(workdir, "config.yaml")

        if not os.path.isfile(config_yaml_path):
            with open(config_yaml_path, "w") as config_yaml:
                yaml.dump(config, config_yaml, default_flow_style=False)

        snakemake_unlock_dir_command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            config_yaml_path,
            "--unlock",
        ]
        run(snakemake_unlock_dir_command, cwd=workdir)

        snakemake_unlock_dir_command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            config_yaml_path,
            "--unlock",
        ]
        run(snakemake_unlock_dir_command, cwd=workdir)

        dag_command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            config_yaml_path,
            "--dag",
            "-n",
        ]
        print_dag = ["dot", "-Tsvg", "-o", "workflow_dag.svg"]

        create_dag = Popen(dag_command, cwd=workdir, stdout=PIPE)
        print_dat_to_file = check_output(
            print_dag, cwd=workdir, stdin=create_dag.stdout
        )
        create_dag.wait()

        snakemake_command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "-j",
            str(AppConfig.THREADS),
            "--configfile",
            config_yaml_path,
        ]

        print(" ".join(snakemake_command))
        run(snakemake_command, cwd=workdir)

        snakemake_report_command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            config_yaml_path,
            "--report",
            "report.html",
        ]
        run(snakemake_report_command, cwd=workdir)
