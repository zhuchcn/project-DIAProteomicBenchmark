"""Parse VEP TSV output to moPepGen GVF."""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import shlex
import subprocess as sp
import sys

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_CODON_TABLE,
    DEFAULT_GVF_OUTPUT,
    DEFAULT_MOPEPGEN_INDEX_DIR,
    DEFAULT_MOPEPGEN_REPO,
    DEFAULT_VEP_OUTPUT_TSV,
)


LOGGER = setup_logger('EcoliCustomDatabase.ParseVEP')


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register ParseVEP arguments."""
    parser.description = 'Parse annotated VEP output for downstream usage.'
    parser.add_argument(
        '--input-tsv',
        type=Path,
        default=DEFAULT_VEP_OUTPUT_TSV,
        help=f'Input VEP TSV path (default: {DEFAULT_VEP_OUTPUT_TSV}).',
    )
    parser.add_argument(
        '--output-gvf',
        type=Path,
        default=DEFAULT_GVF_OUTPUT,
        help=f'Output GVF path (default: {DEFAULT_GVF_OUTPUT}).',
    )
    parser.add_argument(
        '--index-dir',
        type=Path,
        default=DEFAULT_MOPEPGEN_INDEX_DIR,
        help=f'moPepGen index directory (default: {DEFAULT_MOPEPGEN_INDEX_DIR}).',
    )
    parser.add_argument(
        '--codon-table',
        type=str,
        default=DEFAULT_CODON_TABLE,
        help=f'Codon table name for moPepGen parseVEP (default: {DEFAULT_CODON_TABLE}).',
    )
    parser.add_argument(
        '--mopepgen-repo',
        type=Path,
        default=DEFAULT_MOPEPGEN_REPO,
        help=(
            'Path to local moPepGen source repository used via PYTHONPATH '
            f'(default: {DEFAULT_MOPEPGEN_REPO}).'
        ),
    )
    parser.add_argument(
        '--python-exe',
        type=Path,
        default=Path(sys.executable),
        help='Python executable used to run `python -m moPepGen.cli`.',
    )


def run(args: argparse.Namespace) -> None:
    """Execute ParseVEP subcommand."""
    input_tsv = _validate_input_file(args.input_tsv, 'VEP TSV')
    output_gvf = _prepare_output_path(args.output_gvf)
    index_dir = _validate_directory(args.index_dir, 'moPepGen index directory')
    mopepgen_repo = _validate_directory(args.mopepgen_repo, 'moPepGen repository')
    python_exe = _validate_input_file(args.python_exe, 'python executable')

    cmd = [
        str(python_exe),
        '-m',
        'moPepGen.cli',
        'parseVEP',
        '-i',
        str(input_tsv),
        '--source',
        'SNP',
        '-o',
        str(output_gvf),
        '--index-dir',
        str(index_dir),
        '--codon-table',
        args.codon_table,
    ]

    env = os.environ.copy()
    existing_pythonpath = env.get('PYTHONPATH')
    env['PYTHONPATH'] = (
        f'{mopepgen_repo}:{existing_pythonpath}' if existing_pythonpath else str(mopepgen_repo)
    )

    LOGGER.info('Running command: %s', shlex.join(cmd))
    LOGGER.info('Using moPepGen source from: %s', mopepgen_repo)
    sp.run(cmd, check=True, env=env)
    LOGGER.info('Parsed VEP TSV to GVF: %s', output_gvf)


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Validate file existence."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _validate_directory(path_arg: Path, label: str) -> Path:
    """Validate directory existence."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _prepare_output_path(path_arg: Path) -> Path:
    """Create output directory and return resolved output path."""
    path = Path(path_arg).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.resolve()
