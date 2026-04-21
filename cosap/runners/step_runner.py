from concurrent.futures import ProcessPoolExecutor
from functools import partial
from subprocess import CalledProcessError, CompletedProcess, run
from typing import List, Optional

from .._config import AppConfig


def run_command_parallel(
    commands: list, cwd: str = None
) -> List[Optional[CompletedProcess]]:
    max_threads_per_job = AppConfig().MAX_THREADS_PER_JOB

    run_function = partial(run, cwd=cwd, capture_output=True)
    with ProcessPoolExecutor(max_workers=max_threads_per_job) as executor:
        return list(executor.map(run_function, commands))
