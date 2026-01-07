""" Common utilities """
from __future__ import annotations
import os
import logging
from pathlib import Path
from typing import Iterable


MODULE_DIR = Path(__file__).parent
WORK_DIR = MODULE_DIR.parent.parent.parent

## ---- Environ Vars ----
DATA_DIR = os.getenv("DATA_DIR")
assert DATA_DIR is not None, "Please set DATA_DIR in your environment variables."
DATA_DIR = Path(DATA_DIR)

PROJ_DIR = os.getenv("PROJ_DIR")
assert PROJ_DIR is not None, "Please set PROJ_DIR in your environment variables."
PROJ_DIR = Path(PROJ_DIR)

# ---- FragPipe constants --------------------------------------------------
FRAGPIPE_DOCKER_IMAGE = "fcyucn/fragpipe:23.1"
FRAGPIPE_EXE = Path("/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/bin/fragpipe")
PHILOSOPHER_EXE = Path("/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2")

# ---- Files and Directories ----
def get_dataset_dir(dataset_name: str, molecule: str, assay: str) -> Path:
	"""Get the directory path for a specific dataset within the DATA_DIR."""
	dataset_dir = DATA_DIR / 'data' / dataset_name / molecule / assay
	dataset_dir.mkdir(parents=True, exist_ok=True)
	return dataset_dir

def resolve_file(path: Path, label: str) -> Path:
	"""Normalize and ensure a file exists at the provided path."""
	resolved = path.expanduser()
	if not resolved.exists():
		raise FileNotFoundError(f'Missing {label}: {resolved}')
	if not resolved.is_file():
		raise FileNotFoundError(f'Expected {label} to be a file: {resolved}')
	return resolved.resolve()

def resolve_directory(path: Path, label: str, mkdir: bool = False) -> Path:
	"""Normalize and ensure a directory exists at the provided path."""
	resolved = path.expanduser()
	if not resolved.exists():
		if mkdir:
			resolved.mkdir(parents=True, exist_ok=True)
		else:
			raise FileNotFoundError(f'Missing {label}: {resolved}')
	if not resolved.is_dir():
		raise FileNotFoundError(f'Expected {label} to be a directory: {resolved}')
	return resolved.resolve()

def collect_mount_points(*paths: Path) -> list[Path]:
	"""Return ordered directories that should be mounted inside Docker."""
	mounts = {Path.cwd().resolve()}
	for path in paths:
		resolved = path.resolve()
		mounts.add(resolved if resolved.is_dir() else resolved.parent)
	return sorted(mounts)

def build_docker_run_cmd(image: str, mount_points: Iterable[Path], work_dir: Path=None) -> list[str]:
	"""Assemble the docker run prefix for the requested image and mounts."""
	cmd = [
		'docker',
		'run',
		'--rm',
		'-u', f"{os.getuid()}:{os.getgid()}"
	]
	if work_dir is not None:
		cmd.extend(['-w', str(work_dir)])
		if work_dir not in mount_points:
			mount_points = list(mount_points) + [work_dir]
	for mount in mount_points:
		cmd.extend(['-v', f"{mount}:{mount}"])
	cmd.append(image)
	return cmd

# ---- Functions ----
def setup_logger(name=__name__, level=logging.INFO):
	"""Set up and return a logger with a standard format including the command name."""
	logger = logging.getLogger(name)
	if not logger.hasHandlers():
		handler = logging.StreamHandler()
		formatter = logging.Formatter('[%(asctime)s %(name)s] %(message)s')
		handler.setFormatter(formatter)
		logger.addHandler(handler)
	logger.setLevel(level)
	return logger
