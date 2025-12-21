""" Common utilities """
from __future__ import annotations
import os
import logging
from pathlib import Path


MODULE_DIR = Path(__file__).parent
WORK_DIR = MODULE_DIR.parent.parent.parent

## ---- Environ Vars ----
DATA_DIR = os.getenv("DATA_DIR")
assert DATA_DIR is not None, "Please set DATA_DIR in your environment variables."
DATA_DIR = Path(DATA_DIR)

PROJ_DIR = os.getenv("PROJ_DIR")
assert PROJ_DIR is not None, "Please set PROJ_DIR in your environment variables."
PROJ_DIR = Path(PROJ_DIR)

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
