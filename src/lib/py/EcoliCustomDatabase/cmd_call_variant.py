"""Call variant peptides from GVF using moPepGen."""

from __future__ import annotations
import argparse
from pathlib import Path
import shlex

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_CODON_TABLE,
    DEFAULT_GVF_OUTPUT,
    DEFAULT_MOPEPGEN_DOCKER_IMAGE,
    DEFAULT_MOPEPGEN_INDEX_DIR,
    DEFAULT_VARIANT_PEPTIDE_FASTA,
)
from EcoliCustomDatabase.io_utils import run_docker_cmd


LOGGER = setup_logger('EcoliCustomDatabase.CallVariant')


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register CallVariant arguments."""
    parser.description = 'Call moPepGen to generate the custom database FASTA.'
    parser.add_argument(
        '--input-gvf',
        type=Path,
        default=DEFAULT_GVF_OUTPUT,
        help=f'Input GVF file (default: {DEFAULT_GVF_OUTPUT}).',
    )
    parser.add_argument(
        '--output-fasta',
        type=Path,
        default=DEFAULT_VARIANT_PEPTIDE_FASTA,
        help=f'Output FASTA file (default: {DEFAULT_VARIANT_PEPTIDE_FASTA}).',
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
        help=f'Codon table for callVariant (default: {DEFAULT_CODON_TABLE}).',
    )
    parser.add_argument(
        '--docker-image',
        type=str,
        default=DEFAULT_MOPEPGEN_DOCKER_IMAGE,
        help=f'Docker image providing moPepGen CLI (default: {DEFAULT_MOPEPGEN_DOCKER_IMAGE}).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute CallVariant subcommand."""
    input_gvf = _validate_input_file(args.input_gvf, 'input GVF')
    output_fasta = _prepare_output_path(args.output_fasta)
    index_dir = _validate_directory(args.index_dir, 'moPepGen index directory')

    cmd = [
        'moPepGen',
        'callVariant',
        '-i',
        str(input_gvf),
        '-o',
        str(output_fasta),
        '--index-dir',
        str(index_dir),
        '--codon-table',
        args.codon_table,
    ]

    mount_paths = [input_gvf, output_fasta.parent, index_dir]
    LOGGER.info('Running command in Docker: %s', shlex.join(cmd))
    _, docker_cmd = run_docker_cmd(
        image=args.docker_image,
        command=cmd,
        work_dir=output_fasta.parent,
        mount_paths=mount_paths,
    )
    LOGGER.info('Docker command: %s', shlex.join(docker_cmd))
    LOGGER.info('Variant peptide FASTA generated: %s', output_fasta)


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
