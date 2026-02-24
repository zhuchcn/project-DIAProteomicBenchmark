"""Parse VEP TSV output to moPepGen GVF."""

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
    DEFAULT_VEP_OUTPUT_TSV,
)
from EcoliCustomDatabase.io_utils import run_docker_cmd


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
        '--docker-image',
        type=str,
        default=DEFAULT_MOPEPGEN_DOCKER_IMAGE,
        help=f'Docker image providing moPepGen CLI (default: {DEFAULT_MOPEPGEN_DOCKER_IMAGE}).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute ParseVEP subcommand."""
    input_tsv = _validate_input_file(args.input_tsv, 'VEP TSV')
    output_gvf = _prepare_output_path(args.output_gvf)
    index_dir = _validate_directory(args.index_dir, 'moPepGen index directory')

    cmd = [
        'moPepGen',
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

    mount_paths = [input_tsv, output_gvf.parent, index_dir]
    LOGGER.info('Running command in Docker: %s', shlex.join(cmd))
    _, docker_cmd = run_docker_cmd(
        image=args.docker_image,
        command=cmd,
        work_dir=output_gvf.parent,
        mount_paths=mount_paths,
    )
    LOGGER.info('Docker command: %s', shlex.join(docker_cmd))
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
