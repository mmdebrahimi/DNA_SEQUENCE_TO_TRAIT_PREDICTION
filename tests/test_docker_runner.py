"""Unit tests for tools.docker_runner — mocked subprocess, no real Docker."""
from __future__ import annotations

import subprocess
from unittest.mock import patch

import pytest

from tools.docker_runner import DockerRunnerError, run


def _fake_completed(returncode=0, stdout="", stderr=""):
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=stdout, stderr=stderr)


def test_run_builds_docker_command_with_image_and_args():
    with patch("tools.docker_runner.subprocess.run", return_value=_fake_completed()) as m:
        run("mash:2.3", ["dist", "a.fna", "b.fna"])
    cmd = m.call_args[0][0]
    assert cmd[:3] == ["docker", "run", "--rm"]
    assert cmd[3] == "mash:2.3"
    assert cmd[4:] == ["dist", "a.fna", "b.fna"]


def test_run_translates_windows_mounts_to_forward_slash(tmp_path):
    with patch("tools.docker_runner.subprocess.run", return_value=_fake_completed()) as m:
        run(
            "ncbi/amr:4.2.7-2026-03-24.1",
            ["--database", "/db/latest"],
            mounts={str(tmp_path): "/db"},
        )
    cmd = m.call_args[0][0]
    assert "-v" in cmd
    v_idx = cmd.index("-v")
    mount_spec = cmd[v_idx + 1]
    host, container = mount_spec.rsplit(":", 1)
    assert container == "/db"
    assert "\\" not in host


def test_run_passes_env_vars_via_dash_e():
    with patch("tools.docker_runner.subprocess.run", return_value=_fake_completed()) as m:
        run("mash:2.3", ["dist"], env={"FOO": "bar", "BAZ": "qux"})
    cmd = m.call_args[0][0]
    assert "-e" in cmd
    e_specs = [cmd[i + 1] for i, t in enumerate(cmd) if t == "-e"]
    assert "FOO=bar" in e_specs
    assert "BAZ=qux" in e_specs


def test_run_raises_dockerrunnererror_on_nonzero_exit_when_check_true():
    with patch(
        "tools.docker_runner.subprocess.run",
        return_value=_fake_completed(returncode=1, stderr="oops"),
    ):
        with pytest.raises(DockerRunnerError, match="docker exited 1"):
            run("mash:2.3", ["dist"])


def test_run_returns_completed_process_on_nonzero_when_check_false():
    with patch(
        "tools.docker_runner.subprocess.run",
        return_value=_fake_completed(returncode=2, stdout="ok"),
    ):
        proc = run("mash:2.3", ["dist"], check=False)
    assert proc.returncode == 2
    assert proc.stdout == "ok"


def test_run_passes_timeout_to_subprocess():
    with patch("tools.docker_runner.subprocess.run", return_value=_fake_completed()) as m:
        run("mash:2.3", ["dist"], timeout=30.0)
    assert m.call_args.kwargs["timeout"] == 30.0


def test_run_wraps_timeout_expired_as_dockerrunnererror():
    with patch(
        "tools.docker_runner.subprocess.run",
        side_effect=subprocess.TimeoutExpired(cmd=["docker"], timeout=5.0),
    ):
        with pytest.raises(DockerRunnerError, match="exceeded timeout=5.0s"):
            run("mash:2.3", ["dist"], timeout=5.0)


def test_run_wraps_filenotfounderror_as_dockerrunnererror():
    with patch(
        "tools.docker_runner.subprocess.run",
        side_effect=FileNotFoundError(2, "No such file or directory", "docker"),
    ):
        with pytest.raises(DockerRunnerError, match="docker.* not found on PATH"):
            run("mash:2.3", ["dist"])


def test_run_supports_read_only_mount_suffix(tmp_path):
    with patch("tools.docker_runner.subprocess.run", return_value=_fake_completed()) as m:
        run("mash:2.3", ["dist"], mounts={str(tmp_path): "/refseq:ro"})
    cmd = m.call_args[0][0]
    v_idx = cmd.index("-v")
    assert cmd[v_idx + 1].endswith(":/refseq:ro")
