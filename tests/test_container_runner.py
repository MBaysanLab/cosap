import importlib.util
import os
import sys
import types

import pytest


def _load_docker_runner_module():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    module_path = os.path.join(project_root, "cosap", "runners", "docker_runner.py")

    cosap_pkg = types.ModuleType("cosap")
    cosap_pkg.__path__ = [os.path.join(project_root, "cosap")]
    sys.modules.setdefault("cosap", cosap_pkg)

    runners_pkg = types.ModuleType("cosap.runners")
    runners_pkg.__path__ = [os.path.join(project_root, "cosap", "runners")]
    sys.modules.setdefault("cosap.runners", runners_pkg)

    config_module = types.ModuleType("cosap._config")

    class StubAppConfig:
        LIBRARY_PATH = "/refs"
        CONTAINER_RUNTIME = "docker"
        CONTAINER_ENGINE_BIN = "apptainer"

    config_module.AppConfig = StubAppConfig
    sys.modules["cosap._config"] = config_module

    spec = importlib.util.spec_from_file_location(
        "cosap.runners.docker_runner", module_path
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules["cosap.runners.docker_runner"] = module
    spec.loader.exec_module(module)
    return module


docker_runner = _load_docker_runner_module()


class TestContainerRuntimeSelection:
    def test_explicit_docker_runtime(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_RUNTIME", "docker")
        assert docker_runner.get_container_runtime() == "docker"

    def test_explicit_apptainer_runtime(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_RUNTIME", "apptainer")
        assert docker_runner.get_container_runtime() == "apptainer"

    def test_auto_runtime_prefers_apptainer_when_available(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_RUNTIME", "auto")
        monkeypatch.setattr(
            docker_runner,
            "which",
            lambda cmd: "/usr/bin/apptainer" if cmd == "apptainer" else None,
        )
        assert docker_runner.get_container_runtime() == "apptainer"

    def test_auto_runtime_prefers_current_apptainer_container(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_RUNTIME", "auto")
        monkeypatch.setenv("APPTAINER_NAME", "cosap.sif")
        monkeypatch.setattr(docker_runner, "which", lambda cmd: "/usr/bin/docker")

        assert docker_runner.get_container_runtime() == "apptainer"

    def test_invalid_runtime_raises(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_RUNTIME", "podman")
        with pytest.raises(ValueError, match="COSAP_CONTAINER_RUNTIME"):
            docker_runner.get_container_runtime()


class TestApptainerRunner:
    def test_resolve_apptainer_bin_falls_back_to_absolute_candidates(self, monkeypatch):
        monkeypatch.setattr(docker_runner, "which", lambda cmd: None)
        monkeypatch.setattr(docker_runner.os.path, "isfile", lambda path: path == "/usr/bin/apptainer")

        assert docker_runner._resolve_apptainer_bin("apptainer") == "/usr/bin/apptainer"

    def test_inherited_bind_paths_are_preserved(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "LIBRARY_PATH", "/refs")
        monkeypatch.setenv("APPTAINER_BIND", "/outer/work:/outer/work:rw,/shared:/shared:ro")

        runner = docker_runner._ApptainerBackendRunner()
        bind_paths = runner._build_bind_paths(
            "/workdir",
            [{"source": "/input", "target": "/input", "mode": "ro"}],
        )

        assert ("/outer/work", "/outer/work", "rw") in bind_paths
        assert ("/shared", "/shared", "ro") in bind_paths
        assert ("/refs", "/refs", "ro") in bind_paths
        assert ("/workdir", "/workdir", "rw") in bind_paths
        assert ("/input", "/input", "ro") in bind_paths

    def test_apptainer_command_renders_gpu_and_binds(self):
        command = docker_runner._ApptainerBackendRunner._container_run_to_cmd(
            apptainer_bin="apptainer",
            image="google/deepvariant:1.10.0-beta",
            command=["run_deepvariant", "--help"],
            working_dir="/workdir",
            bind_paths=[("/refs", "/refs", "ro"), ("/workdir", "/workdir", "rw")],
            use_gpu=True,
        )

        assert "apptainer exec" in command
        assert "--nv" in command
        assert "--pwd /workdir" in command
        assert "--bind /refs:/refs:ro" in command
        assert "--bind /workdir:/workdir" in command
        assert "docker://google/deepvariant:1.10.0-beta" in command

    def test_apptainer_args_skip_gpu_flag_for_cpu(self):
        command = docker_runner._ApptainerBackendRunner._container_run_to_args(
            apptainer_bin="apptainer",
            image="google/deepvariant:1.10.0-beta",
            command=["run_deepvariant", "--help"],
            working_dir="/workdir",
            bind_paths=[("/refs", "/refs", "ro")],
            use_gpu=False,
        )

        assert command[:2] == ["apptainer", "exec"]
        assert "--nv" not in command
        assert command[-3:] == ["docker://google/deepvariant:1.10.0-beta", "run_deepvariant", "--help"]

    def test_missing_child_apptainer_binary_has_nested_runtime_hint(self, monkeypatch):
        monkeypatch.setenv("APPTAINER_NAME", "cosap.sif")

        runner = docker_runner._ApptainerBackendRunner()
        monkeypatch.setattr(docker_runner, "which", lambda cmd: None)
        monkeypatch.setattr(docker_runner.os.path, "isfile", lambda path: False)

        with pytest.raises(RuntimeError, match="Rebuild the COSAP image"):
            runner._run_impl(
                image="google/deepvariant:1.10.0-beta",
                command=["run_deepvariant", "--help"],
                workdir="/workdir",
                bind_paths=[],
                remove=True,
                run_as_root=False,
            )

    def test_docker_runner_selects_apptainer_backend(self, monkeypatch):
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_RUNTIME", "apptainer")
        monkeypatch.setattr(docker_runner.AppConfig, "CONTAINER_ENGINE_BIN", "apptainer")

        runner = docker_runner.DockerRunner()

        assert runner.runtime == "apptainer"
        assert isinstance(runner._runner, docker_runner._ApptainerBackendRunner)
