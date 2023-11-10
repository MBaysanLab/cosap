import os

import docker

from ..._config import AppConfig


class DockerRunner:
    def __init__(self) -> None:
        self.docker_client = docker.from_env()

    def run(self, image: str, command: str, workdir: str = None) -> None:
        library_path = AppConfig.LIBRARY_PATH

        hostname = os.getenv("HOSTNAME") if self._check_if_running_in_docker() else None

        workdir = (
            self.docker_client.containers.get(hostname).attrs["Config"]["WorkingDir"]
            if self._check_if_running_in_docker()
            else workdir
        )

        volumes = (
            {
                library_path: {"bind": library_path, "mode": "ro"},
                workdir: {"bind": workdir, "mode": "rw"},
            }
            if not self._check_if_running_in_docker()
            else None
        )
        volumes_from = [hostname] if self._check_if_running_in_docker() else None

        self.docker_client.containers.run(
            image=image,
            command=command,
            working_dir=workdir,
            volumes=volumes,
            volumes_from=volumes_from,
            remove=True,
            detach=False,
            restart_policy={"Name": "no"},
            device_requests=[
                docker.types.DeviceRequest(count=-1, capabilities=[["gpu"]])
            ],
        )

    def _check_if_running_in_docker(self) -> bool:
        """
        Returns True if running in docker container.
        """
        return os.path.exists("/.dockerenv")