import multiprocessing
import os
import sys
from subprocess import PIPE, Popen, check_output, run
from typing import Dict, List

import yaml

from .._config import AppConfig
from .._utils import join_paths


class SnakemakeRunner:
    def __init__(self, pipeline_config, workdir):
        self.pipeline_config = pipeline_config
        self.workdir = workdir

    def _create_unlock_dir_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.pipeline_config,
            "--unlock",
        ]
        return command

    def _create_workflow_dag_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.pipeline_config,
            "--dag",
            "-n",
        ]
        return command

    def _create_save_dag_as_svg_command(self) -> list:
        command = ["dot", "-Tsvg", "-o", "workflow_dag.svg"]
        return command

    def _create_dry_run_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.pipeline_config,
            "-r",
            "-n",
        ]
        return command

    def _create_snakemake_run_command(self) -> list:
        available_cpu = multiprocessing.cpu_count()
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "-j",
            str(available_cpu // AppConfig.THREADS),
            "--configfile",
            self.pipeline_config,
            "-r",
            "--use-conda",
        ]
        return command

    def _create_snakemake_report_command(self) -> list:
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "--configfile",
            self.pipeline_config,
            "--report",
            "report.html",
        ]
        return command

    def run_snakemake_pipeline(self):
        unlock_dir = self._create_unlock_dir_command()
        create_dag = self._create_workflow_dag_command()
        save_dag = self._create_save_dag_as_svg_command()
        dry_run = self._create_dry_run_command()
        snakemake = self._create_snakemake_run_command()
        report = self._create_snakemake_report_command()

        run(unlock_dir, cwd=self.workdir)
        dag = Popen(create_dag, cwd=self.workdir, stdout=PIPE)
        print_dat_to_file = check_output(save_dag, cwd=self.workdir, stdin=dag.stdout)
        dag.wait()
        # cont = input(("Check the DAG of the created workflow. Do you want to continue? ([y]/n)") or "y")
        # if cont.lower() == "n":
        #     sys.exit()
        run(dry_run, cwd=self.workdir)
        run(snakemake, cwd=self.workdir)
        # run(report, cwd=self.workdir)
