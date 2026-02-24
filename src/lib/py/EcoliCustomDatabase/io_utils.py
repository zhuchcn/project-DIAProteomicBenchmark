"""Shared helpers for E. coli custom database subcommands."""

from __future__ import annotations
import argparse
import shlex
import shutil
import subprocess as sp
from pathlib import Path

from EcoliCustomDatabase.constants import DEFAULT_CONDA_ENV
from Common import build_docker_run_cmd, collect_mount_points


def add_conda_env_arg(parser: argparse.ArgumentParser) -> None:
    """Add common conda environment argument."""
    parser.add_argument(
        '--conda-env',
        type=str,
        default=DEFAULT_CONDA_ENV,
        help=f'Conda environment used to execute external tools (default: {DEFAULT_CONDA_ENV}).',
    )


def wrap_cmd(cmd: list[str], conda_env: str | None) -> list[str]:
    """Wrap a command with conda run if environment is specified."""
    if conda_env:
        return ['conda', 'run', '--no-capture-output', '-n', conda_env, *cmd]
    return cmd


def format_cmd(cmd: list[str], conda_env: str | None) -> str:
    """Format command for logging."""
    return shlex.join(wrap_cmd(cmd, conda_env))


def run_cmd(cmd: list[str], conda_env: str | None, **kwargs):
    """Run external command, optionally via conda environment."""
    wrapped = wrap_cmd(cmd, conda_env)
    if 'check' not in kwargs:
        kwargs['check'] = True
    return sp.run(wrapped, **kwargs)


def popen_cmd(cmd: list[str], conda_env: str | None, **kwargs):
    """Start external command process, optionally via conda environment."""
    wrapped = wrap_cmd(cmd, conda_env)
    return sp.Popen(wrapped, **kwargs)


def require_tool(tool_name: str, conda_env: str | None) -> None:
    """Validate required runtime tool availability."""
    if conda_env:
        if shutil.which('conda') is None:
            raise FileNotFoundError('Conda executable not found on PATH, but --conda-env was provided.')
        return
    if shutil.which(tool_name) is None:
        raise FileNotFoundError(f'Required tool not found on PATH: {tool_name}')


def run_docker_cmd(
    image: str,
    command: list[str],
    work_dir: Path | None = None,
    mount_paths: list[Path] | None = None,
):
    """Run a command inside Docker with host paths mounted 1:1."""
    if shutil.which('docker') is None:
        raise FileNotFoundError('Docker executable not found on PATH.')
    mounts = collect_mount_points(*(mount_paths or []))
    cmd = build_docker_run_cmd(image=image, mount_points=mounts, work_dir=work_dir)
    cmd.extend(command)
    return sp.run(cmd, check=True), cmd
