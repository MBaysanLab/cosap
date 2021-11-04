import os
from subprocess import PIPE, Popen, check_output, run
from typing import Dict, List

import yaml

from .._config import AppConfig
from .._pipeline_config import PipelineKeys
from .._utils import join_paths


class SnakemakeRunner:
    def __init__(self, pipeline_config):
        self.pipeline_congfig = pipeline_config
        self.workdir = AppConfig.WORKDIR
        self.config_yaml_path = join_paths(self.workdir, "config.yaml")

    def _write_config_to_yaml(self):
        self.config[PipelineKeys.WORKDIR] = self.workdir

        if not os.path.isfile(self.config_yaml_path):
            with open(self.config_yaml_path, "w") as config_yaml:
                yaml.dump(self.config, config_yaml, default_flow_style=False)

    def _create_unlock_dir_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.config_yaml_path,
            "--unlock",
        ]
        return command

    def _create_workflow_dag_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.config_yaml_path,
            "--dag",
            "-n",
        ]
        return command

    def _create_save_dag_as_svg_command(self) -> list:
        command = ["dot", "-Tsvg", "-o", "workflow_dag.svg"]
        return command

    def _create_snakemake_run_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "-j",
            str(AppConfig.THREADS),
            "--configfile",
            self.config_yaml_path,
        ]
        return command

    def _create_snakemake_report_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.config_yaml_path,
            "--report",
            "report.html",
        ]
        return command

    def run_snakemake_pipeline(self):
        unlock_dir = self._create_unlock_dir_command()
        create_dag = self._create_workflow_dag_command()
        save_dag = self._create_save_dag_as_svg_command()
        snakemake = self._create_snakemake_run_command()
        report = self._create_snakemake_report_command()

        run(unlock_dir, cwd=self.workdir)
        dag = Popen(create_dag, cwd=self.workdir, stdout=PIPE)
        print_dat_to_file = check_output(save_dag, cwd=self.workdir, stdin=dag.stdout)
        dag.wait()
        run(snakemake, cwd=self.workdir)
        run(report, cwd=self.workdir)
