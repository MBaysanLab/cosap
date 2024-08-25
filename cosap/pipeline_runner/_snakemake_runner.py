import multiprocessing
from subprocess import PIPE, Popen, check_output, run
import sys
import psutil

from .._config import AppConfig
from ..pipeline_runner.runners._docker_runner import DockerRunner


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
        slurm = AppConfig.SLURM_CLUSTER
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
        if slurm:
            command.append("--slurm")
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

        # Run and return sys output
        try:
            snakemake_process = Popen(snakemake, cwd=self.workdir, text=True, stderr=sys.stderr)
            snakemake_process.wait()
        except KeyboardInterrupt:
            snakemake_process.terminate()

            snakemake_child_processes = psutil.Process(snakemake_process.pid).children(recursive=True)

            for child_process in snakemake_child_processes:
                child_process.kill()  # Child processes do not respond to `terminate`

            snakemake_child_process_pids = [process.pid for process in snakemake_child_processes]

            DockerRunner.stop_all_cosap_handled_containers(pids=snakemake_child_process_pids)

            snakemake_process.wait()
        
        results = snakemake_process

        return {
            "stdout": results.stdout,
            "stderr": results.stderr,
            "returncode": results.returncode,
        }
