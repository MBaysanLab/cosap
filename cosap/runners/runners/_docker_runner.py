import logging
import os
from typing import List, Union

import docker

from ..._config import AppConfig


class DockerRunner:
    def __init__(self, device: str = "cpu") -> None:
        self.device = device
        self.docker_client = docker.from_env()

    def run(
        self,
        image: str,
        command: Union[str, list],
        workdir: str = None,
        paths_to_bind: List[str] = None,
        remove: bool = True,
    ) -> None:
        if paths_to_bind is None:
            paths_to_bind = []

        try:
            # Check if the image exists
            if not self._image_exists(image):
                logging.info(f"Pulling image {image}")
                self._pull_image(image)

            library_path = AppConfig.LIBRARY_PATH
            running_in_docker = self._running_in_docker()
            hostname = os.getenv("HOSTNAME") if running_in_docker else None

            # Set volumes
            volumes_dict = {
                library_path: {"bind": library_path, "mode": "ro"},
                workdir: {"bind": workdir, "mode": "rw"},
            }
            for path in paths_to_bind:
                volumes_dict[path] = {"bind": path, "mode": "rw"}

            volumes = volumes_dict if not running_in_docker else None
            volumes_from = [hostname] if running_in_docker else None

            # Run and return the log generator
            container = self.docker_client.containers.run(
                image=image,
                command=command,
                working_dir=workdir,
                volumes=volumes,
                volumes_from=volumes_from,
                remove=remove,
                detach=True,
                restart_policy={"Name": "no"},
                device_requests=(
                    [docker.types.DeviceRequest(count=-1, capabilities=[["gpu"]])]
                    if self.device == "gpu"
                    else None
                ),
            )
            logs = container.attach(stdout=True, stderr=True, stream=True, logs=True)

            # Print logs
            for log in logs:
                print(log.decode("utf-8"), end="")

        except docker.errors.DockerException as e:
            logging.error(f"Docker error: {e}")
        except Exception as e:
            logging.error(f"Unexpected error: {e}")

    def _running_in_docker(self) -> bool:
        """
        Returns True if running in docker container.
        """
        return os.path.exists("/.dockerenv")

    def _image_exists(self, image: str) -> bool:
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
