import os
from typing import Union

import docker

from ..._config import AppConfig


class DockerRunner:
    container_label_prefix = "cosap_handled_container_"

    def __init__(self, device: str = "cpu") -> None:
        self.device = device
        self.docker_client = docker.from_env()
        self.container_label = DockerRunner.container_label_prefix + str(os.getpid())

    def run(
        self,
        image: str,
        command: Union[str, list],
        workdir: str = None,
        paths_to_bind: list = [],
        remove: bool = True,
    ) -> None:

        # Check if the image exists
        if not self._check_if_image_exists(image):
            self._pull_image(image)

        library_path = AppConfig.LIBRARY_PATH
        hostname = os.getenv("HOSTNAME") if self._check_if_running_in_docker() else None

        # Set volumes
        volumes_dict = {
                library_path: {"bind": library_path, "mode": "ro"},
                workdir: {"bind": workdir, "mode": "rw"},

            }
        for path in paths_to_bind:
            volumes_dict[path] = {"bind": path, "mode": "rw"}

        volumes = None
        volumes_from = None
        device_requests = None

        if not self._check_if_running_in_docker():
            volumes = volumes_dict

        if self._check_if_running_in_docker():
            volumes_from = [hostname]

        if self.device == "gpu":
            device_requests = [docker.types.DeviceRequest(count=-1, capabilities=[["gpu"]])]
        
        # Run and return the log generator
        container = self.docker_client.containers.run(
            image=image,
            command=command,
            working_dir=workdir,
            volumes=volumes,
            volumes_from=volumes_from,
            remove=True,
            detach=True,
            init=True,
            restart_policy={"Name": "no"},
            device_requests=device_requests,
            labels=[self.container_label]
        )

        logs = container.attach(stdout=True, stderr=True, stream=True, logs=True)

        # Print logs
        for log in logs:
            print(log.decode("utf-8"), end="")

    def _check_if_running_in_docker(self) -> bool:
        """
        Returns True if running in docker container.
        """
        return os.path.exists("/.dockerenv")

    def _check_if_image_exists(self, image: str) -> bool:
        """
        Returns True if image exists.
        """
        try:
            self.docker_client.images.get(image)
            return True
        except docker.errors.ImageNotFound:
            return False

    def _pull_image(self, image: str) -> None:
        """
        Pulls the image.
        """
        self.docker_client.images.pull(image)
    
    @classmethod
    def stop_all_cosap_handled_containers(cls, pids: tuple) -> None:
        """
        pids: PIDs of processes which created the Docker containers to be stopped

        Stops all Docker containers that were started by the processes specified

        Intended to be used with a listing of all child processes of a parent process.

        Example usage:
        
        `processes = psutil.Process(parent_process.pid).children(recursive=True)`

        `process_pids = [process.pid for process in processes]`

        `DockerRunner.stop_all_cosap_handled_containers(pids=process_pids)`
        """
        for pid in pids:
            container_label = cls.container_label_prefix + str(pid)

            docker_client = docker.from_env()
            containers = docker_client.containers
            cosap_containers = containers.list(filters={"label": container_label})

            for container in cosap_containers:
                container.stop()
