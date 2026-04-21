import multiprocessing
import os
from subprocess import PIPE, Popen, check_output, run

import psutil

from .._config import AppConfig
from .docker_runner import DockerRunner


class SnakemakeRunner:
    def __init__(self, pipeline_config, workdir, device):
        self.pipeline_config = pipeline_config
        self.workdir = workdir
        self.device = device

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
            "--rerun-incomplete",
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
            "--rerun-incomplete",
        ]
        return command

    def _create_snakemake_run_command(self) -> list:
        available_cpu = multiprocessing.cpu_count()
        # slurm = AppConfig.SLURM_CLUSTER
        command = [
            "snakemake",
            "-s",
            AppConfig.SNAKEFILE_PATH,
            "-j",
            str(available_cpu // AppConfig.MAX_THREADS_PER_JOB),
            "--configfile",
            self.pipeline_config,
            "--config",
            f"device={self.device}",
            "-r",
            "--use-conda",
            "--rerun-incomplete",
        ]
        # if slurm:
        #     command.append("--slurm")
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
        self._create_snakemake_report_command()

        # Create log file path
        log_file_path = os.path.join(self.workdir, "cosap_output.log")

        # Run unlock command and log output
        with open(log_file_path, "w") as log_file:
            log_file.write("=== Snakemake Pipeline Execution Log ===\n\n")
            log_file.write("1. Unlocking directory...\n")
            log_file.flush()

            unlock_result = run(
                unlock_dir, cwd=self.workdir, capture_output=True, text=True
            )
            log_file.write(f"Unlock stdout: {unlock_result.stdout}\n")
            log_file.write(f"Unlock stderr: {unlock_result.stderr}\n\n")

            log_file.write("2. Creating DAG...\n")
            log_file.flush()

        try:
            dag = Popen(create_dag, cwd=self.workdir, stdout=PIPE)
            print_dat_to_file = check_output(
                save_dag, cwd=self.workdir, stdin=dag.stdout
            )
            dag.wait()

            with open(log_file_path, "a") as log_file:
                log_file.write("DAG created and saved as SVG.\n\n")
        except (FileNotFoundError, NotADirectoryError, OSError) as e:
            with open(log_file_path, "a") as log_file:
                log_file.write(f"Warning: Could not generate DAG SVG: {e}\n")
                log_file.write("Install graphviz to enable DAG visualization.\n\n")

        with open(log_file_path, "a") as log_file:
            log_file.write("3. Running dry run...\n")
            log_file.flush()

        dry_run_result = run(dry_run, cwd=self.workdir, capture_output=True, text=True)
        with open(log_file_path, "a") as log_file:
            log_file.write(f"Dry run stdout: {dry_run_result.stdout}\n")
            log_file.write(f"Dry run stderr: {dry_run_result.stderr}\n\n")

            log_file.write("4. Running main Snakemake pipeline...\n")
            log_file.flush()

        # Run and return sys output
        try:
            snakemake_process = Popen(
                snakemake, cwd=self.workdir, text=True, stdout=PIPE, stderr=PIPE
            )
            stdout, stderr = snakemake_process.communicate()

            # Save stdout and stderr to log file
            with open(log_file_path, "a") as log_file:
                log_file.write("Main pipeline stdout:\n")
                log_file.write(stdout)
                log_file.write("\nMain pipeline stderr:\n")
                log_file.write(stderr)
                log_file.write(f"\nReturn code: {snakemake_process.returncode}\n")

        except KeyboardInterrupt:
            snakemake_process.terminate()

            snakemake_child_processes = psutil.Process(snakemake_process.pid).children(
                recursive=True
            )

            for child_process in snakemake_child_processes:
                child_process.kill()  # Child processes do not respond to `terminate`

            snakemake_child_process_pids = [
                process.pid for process in snakemake_child_processes
            ]

            DockerRunner.stop_all_cosap_handled_containers(
                pids=snakemake_child_process_pids
            )

            snakemake_process.wait()

            # Log the interruption
            with open(log_file_path, "a") as log_file:
                log_file.write("\nPipeline execution was interrupted by user.\n")

            stdout, stderr = "", ""

        return {
            "stdout": stdout,
            "stderr": stderr,
            "returncode": snakemake_process.returncode,
            "log_file": log_file_path,
        }
