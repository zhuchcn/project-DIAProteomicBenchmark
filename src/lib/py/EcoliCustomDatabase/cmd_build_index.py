"""Generate moPepGen index using E. coli references."""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import shlex
import subprocess as sp
import sys

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_ANNOTATION_GTF,
    DEFAULT_MOPEPGEN_INDEX_DIR,
    DEFAULT_MOPEPGEN_REPO,
    DEFAULT_PROTEOME_FASTA,
    DEFAULT_REFERENCE_GENOME,
)


LOGGER = setup_logger('EcoliCustomDatabase.BuildIndex')


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register BuildIndex arguments."""
    parser.description = 'Generate moPepGen index using E. coli references.'
    parser.add_argument(
        '--genome',
        type=Path,
        default=DEFAULT_REFERENCE_GENOME,
        help=f'Reference genome FASTA (default: {DEFAULT_REFERENCE_GENOME}).',
    )
    parser.add_argument(
        '--annotation-gtf',
        type=Path,
        default=DEFAULT_ANNOTATION_GTF,
        help=f'Annotation GTF (default: {DEFAULT_ANNOTATION_GTF}).',
    )
    parser.add_argument(
        '--proteome-fasta',
        type=Path,
        default=DEFAULT_PROTEOME_FASTA,
        help=f'Proteome FASTA (default: {DEFAULT_PROTEOME_FASTA}).',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=DEFAULT_MOPEPGEN_INDEX_DIR,
        help=f'Output directory for moPepGen index (default: {DEFAULT_MOPEPGEN_INDEX_DIR}).',
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
    parser.add_argument(
        '--codon-table',
        type=str,
        default='Bacterial',
        help='Codon table name for translation in moPepGen (default: Bacterial).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute BuildIndex subcommand."""
    genome = _validate_input_file(args.genome, 'reference genome FASTA')
    annotation_gtf = _validate_input_file(args.annotation_gtf, 'annotation GTF')
    proteome_fasta = _validate_input_file(args.proteome_fasta, 'proteome FASTA')
    output_dir = _ensure_output_dir(args.output_dir)
    mopepgen_repo = _validate_directory(args.mopepgen_repo, 'moPepGen repository')
    python_exe = _validate_input_file(args.python_exe, 'python executable')

    cmd = [
        str(python_exe),
        '-m',
        'moPepGen.cli',
        'generateIndex',
        '-g',
        str(genome),
        '-a',
        str(annotation_gtf),
        '-p',
        str(proteome_fasta),
        '-o',
        str(output_dir),
        '--codon-table',
        args.codon_table,
        '--force'
    ]

    env = os.environ.copy()
    existing_pythonpath = env.get('PYTHONPATH')
    env['PYTHONPATH'] = (
        f'{mopepgen_repo}:{existing_pythonpath}' if existing_pythonpath else str(mopepgen_repo)
    )

    LOGGER.info('Running command: %s', shlex.join(cmd))
    LOGGER.info('Using moPepGen source from: %s', mopepgen_repo)
    sp.run(cmd, check=True, env=env)
    LOGGER.info('moPepGen index generation completed: %s', output_dir)


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


def _ensure_output_dir(path_arg: Path) -> Path:
    """Create output directory if needed."""
    path = Path(path_arg).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    return path.resolve()
