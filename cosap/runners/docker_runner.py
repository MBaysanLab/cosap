import os
import re
import shlex
import socket
import subprocess
import time
import uuid
from abc import ABC, abstractmethod
from logging import getLogger
from shutil import which
from typing import Dict, List, Optional, Tuple, Union

from .._config import AppConfig

logger = getLogger(__name__)

# Shell metacharacters that should not appear in paths bound to containers.
_DANGEROUS_PATH_CHARS = re.compile(r"[;\|`$\(\)\{\}!&<>]")
_APPTAINER_ENV_VARS = (
    "APPTAINER_NAME",
    "APPTAINER_CONTAINER",
    "SINGULARITY_NAME",
    "SINGULARITY_CONTAINER",
)
_DEFAULT_APPTAINER_BIN_CANDIDATES = (
    "/usr/bin/apptainer",
    "/usr/local/bin/apptainer",
    "/opt/apptainer/bin/apptainer",
)


def _check_if_running_in_apptainer() -> bool:
    return any(os.getenv(env_var) for env_var in _APPTAINER_ENV_VARS)


def _resolve_apptainer_bin(apptainer_bin: str) -> Optional[str]:
    resolved_bin = which(apptainer_bin)
    if resolved_bin is not None:
        return resolved_bin

    if os.path.isabs(apptainer_bin) and os.path.isfile(apptainer_bin):
        return apptainer_bin

    for candidate in _DEFAULT_APPTAINER_BIN_CANDIDATES:
        if os.path.isfile(candidate):
            return candidate

    return None


class _ContainerRunnerBase(ABC):
    container_label_prefix = "cosap_handled_container_"

    def __init__(self, device: str = "cpu") -> None:
        self.device = device
        self.container_label = self.container_label_prefix + str(os.getpid())

    def run(
        self,
        image: str,
        command: Union[str, list],
        workdir: str = None,
        paths_to_bind: List[Dict[str, str]] = None,
        remove: bool = True,
        run_as_root: bool = False,
    ) -> int:
        """
        image: Container image to use
        command: Command to run in the container
        workdir: Working directory to use in the container. If None, uses current working
                    directory.
        paths_to_bind: List of dictionaries specifying additional paths to bind mount.
                        Each dictionary should have 'source' and 'target' keys, and an optional
                        'mode' key (default is 'rw').
        remove: Whether to remove the container after it exits (default is True)
        """
        if paths_to_bind is None:
            paths_to_bind = []

        validated_workdir = self._validate_workdir(workdir)
        bind_paths = self._build_bind_paths(validated_workdir, paths_to_bind)
        return self._run_impl(
            image=image,
            command=command,
            workdir=validated_workdir,
            bind_paths=bind_paths,
            remove=remove,
            run_as_root=run_as_root,
        )

    @classmethod
    @abstractmethod
    def stop_all_cosap_handled_containers(cls, pids: tuple) -> None:
        raise NotImplementedError

    def _validate_workdir(self, workdir: Optional[str]) -> Optional[str]:
        if not workdir:
            return workdir

        workdir = os.path.realpath(workdir)
        if _DANGEROUS_PATH_CHARS.search(workdir):
            raise ValueError(f"Container workdir contains invalid characters: {workdir}")
        return workdir

    def _build_bind_paths(
        self, workdir: Optional[str], paths_to_bind: List[Dict[str, str]]
    ) -> List[Tuple[str, str, str]]:
        bind_paths = []
        library_path = os.path.realpath(AppConfig.LIBRARY_PATH)
        bind_paths.append((library_path, library_path, "ro"))

        if workdir:
            bind_paths.append((workdir, workdir, "rw"))

        for path_binding in paths_to_bind:
            source = os.path.realpath(path_binding["source"])
            target = path_binding["target"]
            mode = path_binding.get("mode", "rw")

            if _DANGEROUS_PATH_CHARS.search(source):
                raise ValueError(
                    f"Container bind mount source contains invalid characters: {source}"
                )

            bind_paths.append((source, target, mode))

        return self._dedupe_bind_paths(bind_paths)

    @staticmethod
    def _dedupe_bind_paths(
        bind_paths: List[Tuple[str, str, str]]
    ) -> List[Tuple[str, str, str]]:
        deduped = []
        seen = set()
        for bind_path in bind_paths:
            if bind_path in seen:
                continue
            seen.add(bind_path)
            deduped.append(bind_path)
        return deduped

    @abstractmethod
    def _run_impl(
        self,
        image: str,
        command: Union[str, list],
        workdir: Optional[str],
        bind_paths: List[Tuple[str, str, str]],
        remove: bool,
        run_as_root: bool,
    ) -> int:
        raise NotImplementedError


class _DockerBackendRunner(_ContainerRunnerBase):
    def __init__(self, device: str = "cpu") -> None:
        super().__init__(device=device)
        self._docker = self._load_docker_module()
        self.docker_client = self._docker.from_env()

    @staticmethod
    def _load_docker_module():
        try:
            import docker
        except ImportError as exc:
            raise RuntimeError(
                "Docker backend selected, but the 'docker' Python package is not installed."
            ) from exc
        return docker

    def _run_impl(
        self,
        image: str,
        command: Union[str, list],
        workdir: Optional[str],
        bind_paths: List[Tuple[str, str, str]],
        remove: bool,
        run_as_root: bool,
    ) -> int:
        if not self._check_if_image_exists(image):
            self._pull_image(image)

        hostname = os.getenv("HOSTNAME") if self._check_if_running_in_docker() else None

        volumes = None
        volumes_from = None
        device_requests = None

        volumes_dict = {
            source: {"bind": target, "mode": mode}
            for source, target, mode in bind_paths
        }

        if not self._check_if_running_in_docker():
            volumes = volumes_dict
        else:
            volumes_from = [hostname]

        if self.device == "gpu":
            device_requests = [
                self._docker.types.DeviceRequest(count=-1, capabilities=[["gpu"]])
            ]

        user = None if run_as_root else os.getuid()

        if run_as_root:
            logger.warning(
                "Running Docker container as root. "
                "This is a security risk when combined with host bind mounts. "
                "Consider running as a non-root user instead."
            )

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
        print(f"Here is the full docker run command:\n{docker_run_command}")

        container = self.docker_client.containers.run(
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

        try:
            logs = container.attach(stdout=True, stderr=True, stream=True, logs=True)

            for log in logs:
                print(log.decode("utf-8").strip())

            wait_result = container.wait()
            status_code = int(wait_result.get("StatusCode", 1))
            return status_code
        except Exception:
            try:
                container.stop(timeout=10)
            except Exception:
                pass
            try:
                container.remove(force=True)
            except Exception:
                pass
            raise

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
        except self._docker.errors.ImageNotFound:
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
        """
        docker = cls._load_docker_module()
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

        if detach:
            cmd_parts.append("-d")

        if init:
            cmd_parts.append("--init")

        if remove:
            cmd_parts.append("--rm")

        if working_dir:
            cmd_parts.extend(["-w", working_dir])

        if user is not None:
            cmd_parts.extend(["-u", str(user)])

        if volumes:
            for host_path, container_info in volumes.items():
                bind = container_info["bind"]
                mode = container_info.get("mode", "rw")
                cmd_parts.extend(["-v", f"{host_path}:{bind}:{mode}"])

        if volumes_from:
            if isinstance(volumes_from, list):
                for vol in volumes_from:
                    cmd_parts.extend(["--volumes-from", vol])
            else:
                cmd_parts.extend(["--volumes-from", volumes_from])

        if restart_policy:
            policy_name = restart_policy.get("Name", "no")
            cmd_parts.extend(["--restart", policy_name])

        if device_requests:
            for device_req in device_requests:
                driver = getattr(device_req, "driver", None)
                if driver is None and hasattr(device_req, "get"):
                    driver = device_req.get("driver")
                capabilities = []
                if hasattr(device_req, "get"):
                    capabilities = device_req.get("capabilities", [])
                else:
                    capabilities = getattr(device_req, "capabilities", [])

                if driver == "nvidia":
                    if capabilities:
                        caps = ",".join(capabilities[0]) if capabilities else "all"
                        cmd_parts.extend(["--gpus", f"all,capabilities={caps}"])
                    else:
                        cmd_parts.extend(["--gpus", "all"])
                elif capabilities:
                    caps = ",".join(capabilities[0]) if capabilities else "gpu"
                    cmd_parts.extend(["--gpus", f"all,capabilities={caps}"])

        if labels:
            if isinstance(labels, list):
                for label in labels:
                    cmd_parts.extend(["--label", label])
            elif isinstance(labels, dict):
                for key, value in labels.items():
                    cmd_parts.extend(["--label", f"{key}={value}"])

        cmd_parts.append(image)

        if command:
            if isinstance(command, list):
                cmd_parts.extend(command)
            else:
                cmd_parts.append(command)

        return " ".join(shlex.quote(part) for part in cmd_parts)


class _ApptainerBackendRunner(_ContainerRunnerBase):
    def __init__(self, device: str = "cpu") -> None:
        super().__init__(device=device)
        self.apptainer_bin = AppConfig.CONTAINER_ENGINE_BIN

    def _run_impl(
        self,
        image: str,
        command: Union[str, list],
        workdir: Optional[str],
        bind_paths: List[Tuple[str, str, str]],
        remove: bool,
        run_as_root: bool,
    ) -> int:
        resolved_apptainer_bin = _resolve_apptainer_bin(self.apptainer_bin)
        if resolved_apptainer_bin is None:
            if _check_if_running_in_apptainer():
                raise RuntimeError(
                    "Apptainer backend selected inside an Apptainer container, but the "
                    f"child Apptainer executable was not found on PATH: {self.apptainer_bin}. "
                    "Rebuild the COSAP image so its bundled Apptainer install is present, or bind a full "
                    "host Apptainer installation into the COSAP container and set COSAP_APPTAINER_BIN to "
                    "the mounted executable. Mounting only the binary is not sufficient because Apptainer "
                    "also needs its config and support files."
                )
            raise RuntimeError(
                f"Apptainer backend selected, but executable was not found: {self.apptainer_bin}"
            )
        self.apptainer_bin = resolved_apptainer_bin

        if run_as_root:
            logger.warning(
                "Apptainer backend does not support running child containers as root. "
                "Continuing as the invoking user."
            )

        if remove is False:
            logger.debug(
                "Apptainer backend ignores the remove flag because exec launches a foreground process."
            )

        apptainer_command = self._container_run_to_cmd(
            apptainer_bin=self.apptainer_bin,
            image=image,
            command=command,
            working_dir=workdir,
            bind_paths=bind_paths,
            use_gpu=self.device == "gpu",
        )
        print(f"Here is the full apptainer exec command:\n{apptainer_command}")

        process_command = self._container_run_to_args(
            apptainer_bin=self.apptainer_bin,
            image=image,
            command=command,
            working_dir=workdir,
            bind_paths=bind_paths,
            use_gpu=self.device == "gpu",
        )

        process = subprocess.Popen(
            process_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        try:
            assert process.stdout is not None
            for line in process.stdout:
                print(line.rstrip())
            return process.wait()
        except Exception:
            try:
                process.terminate()
            except Exception:
                pass
            raise

    def _build_bind_paths(
        self, workdir: Optional[str], paths_to_bind: List[Dict[str, str]]
    ) -> List[Tuple[str, str, str]]:
        bind_paths = super()._build_bind_paths(workdir, paths_to_bind)
        inherited_bind_paths = self._parse_inherited_bind_paths(
            os.getenv("APPTAINER_BIND") or os.getenv("APPTAINER_BINDPATH")
        )
        return self._dedupe_bind_paths(inherited_bind_paths + bind_paths)

    @staticmethod
    def _parse_inherited_bind_paths(bind_env: Optional[str]) -> List[Tuple[str, str, str]]:
        if not bind_env:
            return []

        bind_paths = []
        for raw_spec in bind_env.split(","):
            raw_spec = raw_spec.strip()
            if not raw_spec:
                continue

            parts = raw_spec.split(":")
            source = parts[0]
            target = parts[1] if len(parts) >= 2 and parts[1] else parts[0]
            mode = parts[2] if len(parts) >= 3 and parts[2] else "rw"
            bind_paths.append((source, target, mode))

        return bind_paths

    @staticmethod
    def _image_to_apptainer_uri(image: str) -> str:
        if image.startswith(("docker://", "oras://", "library://")):
            return image
        return f"docker://{image}"

    @classmethod
    def _container_run_to_args(
        cls,
        apptainer_bin: str,
        image: str,
        command=None,
        working_dir=None,
        bind_paths=None,
        use_gpu: bool = False,
    ) -> List[str]:
        cmd_parts = [apptainer_bin, "exec"]

        if use_gpu:
            cmd_parts.append("--nv")

        if working_dir:
            cmd_parts.extend(["--pwd", working_dir])

        if bind_paths:
            for source, target, mode in bind_paths:
                bind_spec = f"{source}:{target}"
                if mode == "ro":
                    bind_spec = f"{bind_spec}:ro"
                cmd_parts.extend(["--bind", bind_spec])

        cmd_parts.append(cls._image_to_apptainer_uri(image))

        if command:
            if isinstance(command, list):
                cmd_parts.extend(command)
            else:
                cmd_parts.append(command)

        return cmd_parts

    @classmethod
    def _container_run_to_cmd(
        cls,
        apptainer_bin: str,
        image: str,
        command=None,
        working_dir=None,
        bind_paths=None,
        use_gpu: bool = False,
    ) -> str:
        return " ".join(
            shlex.quote(part)
            for part in cls._container_run_to_args(
                apptainer_bin=apptainer_bin,
                image=image,
                command=command,
                working_dir=working_dir,
                bind_paths=bind_paths,
                use_gpu=use_gpu,
            )
        )

    @classmethod
    def stop_all_cosap_handled_containers(cls, pids: tuple) -> None:
        # Nested Apptainer executions are foreground child processes, so normal
        # child-process termination is sufficient and there is no daemon-managed
        # container state to clean up here.
        return None


class _ExecutorBackendRunner(_ContainerRunnerBase):
    def __init__(self, device: str = "cpu") -> None:
        super().__init__(device=device)
        try:
            import requests
        except ImportError as exc:
            raise RuntimeError(
                "Executor backend selected, but the 'requests' package is not installed."
            ) from exc

        if not AppConfig.EXECUTOR_TOKEN:
            raise RuntimeError("COSAP_EXECUTOR_TOKEN must be set for executor mode.")

        self._requests = requests
        self._session = requests.Session()
        self._base_url = AppConfig.EXECUTOR_URL
        self._headers = {
            "Authorization": f"Bearer {AppConfig.EXECUTOR_TOKEN}",
            "Content-Type": "application/json",
        }
        self._timeout = (
            AppConfig.EXECUTOR_CONNECT_TIMEOUT,
            AppConfig.EXECUTOR_READ_TIMEOUT,
        )

    def _request(self, method: str, path: str, *, deadline: float = None, **kwargs):
        retry_deadline = time.monotonic() + AppConfig.EXECUTOR_MAX_WAIT_SECONDS
        if deadline is not None:
            retry_deadline = min(retry_deadline, deadline)
        request_timeout = kwargs.pop("timeout", self._timeout)
        while True:
            try:
                response = self._session.request(
                    method,
                    f"{self._base_url}{path}",
                    headers=self._headers,
                    timeout=request_timeout,
                    **kwargs,
                )
            except self._requests.RequestException as exc:
                if time.monotonic() >= retry_deadline:
                    raise RuntimeError(f"COSAP executor request failed: {exc}") from exc
                time.sleep(
                    min(
                        AppConfig.EXECUTOR_RETRY_SECONDS,
                        max(0, retry_deadline - time.monotonic()),
                    )
                )
                continue
            retryable = response.status_code == 429 or response.status_code >= 500
            if not retryable or time.monotonic() >= retry_deadline:
                break
            time.sleep(
                min(
                    AppConfig.EXECUTOR_RETRY_SECONDS,
                    max(0, retry_deadline - time.monotonic()),
                )
            )

        if not response.ok:
            try:
                detail = response.json().get("detail", response.text)
            except ValueError:
                detail = response.text
            raise RuntimeError(
                f"COSAP executor rejected request ({response.status_code}): {detail}"
            )
        return response

    def _run_impl(
        self,
        image: str,
        command: Union[str, list],
        workdir: Optional[str],
        bind_paths: List[Tuple[str, str, str]],
        remove: bool,
        run_as_root: bool,
    ) -> int:
        if run_as_root:
            raise RuntimeError("Executor policy does not permit root tool containers.")

        execution_id = str(uuid.uuid4())
        payload = {
            "execution_id": execution_id,
            "owner": {"hostname": socket.gethostname(), "pid": os.getpid()},
            "image": image,
            "command": command,
            "workdir": workdir,
            "mounts": [
                {"source": source, "target": target, "mode": mode}
                for source, target, mode in bind_paths
            ],
            "device": self.device,
            "lease_seconds": AppConfig.EXECUTOR_LEASE_SECONDS,
        }
        self._request("POST", "/v1/executions", json=payload)

        cursor = 0
        lease_interval = max(1, AppConfig.EXECUTOR_LEASE_SECONDS // 3)
        next_lease = time.monotonic() + lease_interval
        execution_deadline = time.monotonic() + AppConfig.EXECUTOR_MAX_RUN_SECONDS
        terminal = False
        try:
            while True:
                if time.monotonic() >= execution_deadline:
                    raise RuntimeError(
                        f"COSAP executor job exceeded {AppConfig.EXECUTOR_MAX_RUN_SECONDS:g} seconds"
                    )
                if time.monotonic() >= next_lease:
                    lease_deadline = min(
                        execution_deadline,
                        time.monotonic() + max(1, AppConfig.EXECUTOR_LEASE_SECONDS // 2),
                    )
                    self._request(
                        "POST",
                        f"/v1/executions/{execution_id}/lease",
                        json={"lease_seconds": AppConfig.EXECUTOR_LEASE_SECONDS},
                        deadline=lease_deadline,
                    )
                    next_lease = time.monotonic() + lease_interval

                poll_deadline = min(execution_deadline, next_lease)
                try:
                    logs = self._request(
                        "GET",
                        f"/v1/executions/{execution_id}/logs",
                        params={"cursor": cursor},
                        deadline=poll_deadline,
                    ).json()
                    for chunk in logs.get("chunks", logs.get("lines", [])):
                        print(chunk, end="" if "chunks" in logs else "\n")
                    cursor = logs.get("next_cursor", cursor)

                    status = self._request(
                        "GET",
                        f"/v1/executions/{execution_id}",
                        deadline=poll_deadline,
                    ).json()
                except RuntimeError:
                    logger.warning(
                        "Transient executor polling failure for job %s",
                        execution_id,
                        exc_info=True,
                    )
                    time.sleep(
                        min(
                            AppConfig.EXECUTOR_POLL_SECONDS,
                            max(0, next_lease - time.monotonic()),
                        )
                    )
                    continue
                if status["status"] in {"succeeded", "failed", "cancelled"}:
                    while True:
                        if time.monotonic() >= execution_deadline:
                            raise RuntimeError(
                                "COSAP executor timed out while draining final logs"
                            )
                        final_logs = self._request(
                            "GET",
                            f"/v1/executions/{execution_id}/logs",
                            params={"cursor": cursor},
                            deadline=execution_deadline,
                        ).json()
                        chunks = final_logs.get(
                            "chunks", final_logs.get("lines", [])
                        )
                        for chunk in chunks:
                            print(
                                chunk,
                                end="" if "chunks" in final_logs else "\n",
                            )
                        next_cursor = final_logs.get("next_cursor", cursor)
                        if not chunks or next_cursor == cursor:
                            break
                        cursor = next_cursor
                    terminal = True
                    exit_code = status.get("exit_code")
                    if exit_code is None:
                        exit_code = 0 if status["status"] == "succeeded" else 1
                    return int(exit_code)

                time.sleep(AppConfig.EXECUTOR_POLL_SECONDS)
        finally:
            if not terminal:
                try:
                    cancel_timeout = 5
                    self._request(
                        "DELETE",
                        f"/v1/executions/{execution_id}",
                        deadline=time.monotonic() + cancel_timeout,
                        timeout=(
                            min(AppConfig.EXECUTOR_CONNECT_TIMEOUT, cancel_timeout),
                            min(AppConfig.EXECUTOR_READ_TIMEOUT, cancel_timeout),
                        ),
                    )
                except RuntimeError:
                    logger.exception("Failed to cancel executor job %s", execution_id)

    @classmethod
    def stop_all_cosap_handled_containers(cls, pids: tuple) -> None:
        runner = cls()
        runner._request(
            "POST",
            "/v1/executions/cancel-by-owner",
            json={"hostname": socket.gethostname(), "pids": list(pids)},
        )


def get_container_runtime() -> str:
    runtime = AppConfig.CONTAINER_RUNTIME.lower()
    if runtime == "auto":
        if _check_if_running_in_apptainer():
            return "apptainer"
        if os.path.exists("/.dockerenv") or which("docker") is not None:
            return "docker"
        if which(AppConfig.CONTAINER_ENGINE_BIN) is not None:
            return "apptainer"
        return "docker"

    if runtime not in {"docker", "apptainer", "executor"}:
        raise ValueError(
            "COSAP_CONTAINER_RUNTIME must be one of: docker, apptainer, executor, auto"
        )
    return runtime


class DockerRunner:
    container_label_prefix = _ContainerRunnerBase.container_label_prefix

    def __init__(self, device: str = "cpu") -> None:
        self.device = device
        self.runtime = get_container_runtime()
        if self.runtime == "docker":
            self._runner = _DockerBackendRunner(device=device)
        elif self.runtime == "apptainer":
            self._runner = _ApptainerBackendRunner(device=device)
        elif self.runtime == "executor":
            self._runner = _ExecutorBackendRunner(device=device)
        else:
            raise ValueError(f"Unsupported container runtime: {self.runtime}")

    def run(
        self,
        image: str,
        command: Union[str, list],
        workdir: str = None,
        paths_to_bind: List[Dict[str, str]] = None,
        remove: bool = True,
        run_as_root: bool = False,
    ) -> int:
        return self._runner.run(
            image=image,
            command=command,
            workdir=workdir,
            paths_to_bind=paths_to_bind,
            remove=remove,
            run_as_root=run_as_root,
        )

    @classmethod
    def stop_all_cosap_handled_containers(cls, pids: tuple) -> None:
        runtime = get_container_runtime()
        if runtime == "docker":
            _DockerBackendRunner.stop_all_cosap_handled_containers(pids=pids)
            return
        if runtime == "executor":
            _ExecutorBackendRunner.stop_all_cosap_handled_containers(pids=pids)
            return
        _ApptainerBackendRunner.stop_all_cosap_handled_containers(pids=pids)

    @classmethod
    def _container_run_to_cmd(cls, *args, **kwargs):
        return _DockerBackendRunner._container_run_to_cmd(*args, **kwargs)
