"""Install VEP cache for E. coli."""

from __future__ import annotations
import argparse
from pathlib import Path

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_VEP_ASSEMBLY,
    DEFAULT_VEP_CACHE_DIR,
    DEFAULT_VEP_KINGDOM,
    DEFAULT_VEP_SPECIES,
)
from EcoliCustomDatabase.io_utils import add_conda_env_arg, format_cmd, require_tool, run_cmd


LOGGER = setup_logger('EcoliCustomDatabase.InstallVepCache')


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register InstallVepCache arguments."""
    parser.description = 'Install VEP cache for E. coli.'
    add_conda_env_arg(parser)
    parser.add_argument(
        '--species',
        type=str,
        default=DEFAULT_VEP_SPECIES,
        help=f'VEP species argument (default: {DEFAULT_VEP_SPECIES}).',
    )
    parser.add_argument(
        '--kingdom',
        type=str,
        default=DEFAULT_VEP_KINGDOM,
        help=f'VEP kingdom argument (default: {DEFAULT_VEP_KINGDOM}).',
    )
    parser.add_argument(
        '--cache-dir',
        type=Path,
        default=DEFAULT_VEP_CACHE_DIR,
        help=f'VEP cache directory (default: {DEFAULT_VEP_CACHE_DIR}).',
    )
    parser.add_argument(
        '--assembly',
        type=str,
        default=DEFAULT_VEP_ASSEMBLY,
        help=f'Assembly name (default: {DEFAULT_VEP_ASSEMBLY}).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute InstallVepCache subcommand."""
    require_tool('vep_install', args.conda_env)

    cache_dir = Path(args.cache_dir).expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        'vep_install',
        '-a',
        'cf',
        '-s',
        str(args.species),
        '-y',
        str(args.kingdom),
        '-c',
        str(cache_dir),
        '--assembly',
        str(args.assembly),
    ]
    LOGGER.info('Running command: %s', format_cmd(cmd, args.conda_env))
    run_cmd(cmd, args.conda_env)
    LOGGER.info('VEP cache install completed under %s', cache_dir)
