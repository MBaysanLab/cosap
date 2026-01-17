import os
from logging import getLogger
from typing import Dict, List, Union

import docker

from .._config import AppConfig

logger = getLogger(__name__)


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
        paths_to_bind: List[Dict[str, str]] = None,
        remove: bool = True,
        run_as_root: bool = False,
    ) -> None:
        """
        image: Docker image to use
        command: Command to run in the Docker container
        workdir: Working directory to use in the Docker container. If None, uses current working
                    directory.
        paths_to_bind: List of dictionaries specifying additional paths to bind mount.
                        Each dictionary should have 'source' and 'target' keys, and an optional
                        'mode' key (default is 'rw').
        remove: Whether to remove the container after it exits (default is True)
        """

        if paths_to_bind is None:
            paths_to_bind = []

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
        for path_binding in paths_to_bind:
            source = path_binding["source"]
            target = path_binding["target"]
            mode = path_binding.get("mode", "rw")
            volumes_dict[source] = {"bind": target, "mode": mode}

        volumes = None
        volumes_from = None
        device_requests = None

        if not self._check_if_running_in_docker():
            volumes = volumes_dict
        else:
            volumes_from = [hostname]

        if self.device == "gpu":
            device_requests = [
                docker.types.DeviceRequest(count=-1, capabilities=[["gpu"]])
            ]

        user = None if run_as_root else os.getuid()
        # Run and return the log generator
        docker_run_command = self._container_run_to_cmd(
            image=image,
            command=command,
            working_dir=workdir,
            volumes=volumes,
            volumes_from=volumes_from,
            detach=True,
            init=True,
            restart_policy={"Name": "no"},
            device_requests=device_requests,
            labels=[self.container_label],
            remove=remove,
            user=user,
        )
        print(
            f"Here is the full docker run command:\n{docker_run_command}"
        )

        container = self.docker_client.containers.run(
            image=image,
            command=command,
            working_dir=workdir,
            volumes=volumes,
            volumes_from=volumes_from,
            detach=True,
            init=True,  # Required to handle process signals correctly
            restart_policy={"Name": "no"},
            device_requests=device_requests,
            labels=[self.container_label],
            remove=remove,
            user=user,
        )

        logs = container.attach(stdout=True, stderr=True, stream=True, logs=True)

        for log in logs:
            print(log.decode("utf-8").strip())

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

    @classmethod
    def _container_run_to_cmd(
        cls,
        image,
        command=None,
        working_dir=None,
        volumes=None,
        volumes_from=None,
        detach=True,
        init=True,
        restart_policy=None,
        device_requests=None,
        labels=None,
        remove=False,
        user=None,
    ):
        """Convert container.run() parameters to docker run command"""
        cmd_parts = ["docker", "run"]

        # Detach mode
        if detach:
            cmd_parts.append("-d")

        # Init process
        if init:
            cmd_parts.append("--init")

        # Remove container after exit
        if remove:
            cmd_parts.append("--rm")

        # Working directory
        if working_dir:
            cmd_parts.extend(["-w", working_dir])

        # User
        if user is not None:
            cmd_parts.extend(["-u", str(user)])

        # Volumes
        if volumes:
            for host_path, container_info in volumes.items():
                bind = container_info["bind"]
                mode = container_info.get("mode", "rw")
                cmd_parts.extend(["-v", f"{host_path}:{bind}:{mode}"])

        # Volumes from
        if volumes_from:
            if isinstance(volumes_from, list):
                for vol in volumes_from:
                    cmd_parts.extend(["--volumes-from", vol])
            else:
                cmd_parts.extend(["--volumes-from", volumes_from])

        # Restart policy
        if restart_policy:
            policy_name = restart_policy.get("Name", "no")
            cmd_parts.extend(["--restart", policy_name])

        # Device requests (GPU support)
        if device_requests:
            for device_req in device_requests:
                if device_req.get("driver") == "nvidia":
                    # Use --gpus flag for NVIDIA GPUs
                    capabilities = device_req.get("capabilities", [])
                    if capabilities:
                        caps = ",".join(capabilities[0]) if capabilities else "all"
                        cmd_parts.extend(["--gpus", f"all,capabilities={caps}"])
                    else:
                        cmd_parts.extend(["--gpus", "all"])

        # Labels
        if labels:
            if isinstance(labels, list):
                for label in labels:
                    cmd_parts.extend(["--label", label])
            elif isinstance(labels, dict):
                for key, value in labels.items():
                    cmd_parts.extend(["--label", f"{key}={value}"])

        # Image
        cmd_parts.append(image)

        # Command
        if command:
            if isinstance(command, list):
                cmd_parts.extend(command)
            else:
                cmd_parts.append(command)

        return " ".join(cmd_parts)
