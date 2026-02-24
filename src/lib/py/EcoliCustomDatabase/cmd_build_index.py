"""Generate moPepGen index using E. coli references."""

from __future__ import annotations
import argparse
from pathlib import Path
import shlex

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_ANNOTATION_GTF,
    DEFAULT_MOPEPGEN_DOCKER_IMAGE,
    DEFAULT_MOPEPGEN_INDEX_DIR,
    DEFAULT_PROTEOME_FASTA,
    DEFAULT_REFERENCE_GENOME,
)
from EcoliCustomDatabase.io_utils import run_docker_cmd


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
        '--docker-image',
        type=str,
        default=DEFAULT_MOPEPGEN_DOCKER_IMAGE,
        help=f'Docker image providing moPepGen CLI (default: {DEFAULT_MOPEPGEN_DOCKER_IMAGE}).',
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

    cmd = [
        'moPepGen',
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

    mount_paths = [genome, annotation_gtf, proteome_fasta, output_dir]
    LOGGER.info('Running command in Docker: %s', shlex.join(cmd))
    _, docker_cmd = run_docker_cmd(
        image=args.docker_image,
        command=cmd,
        work_dir=output_dir,
        mount_paths=mount_paths,
    )
    LOGGER.info('Docker command: %s', shlex.join(docker_cmd))
    LOGGER.info('moPepGen index generation completed: %s', output_dir)


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Validate file existence."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _ensure_output_dir(path_arg: Path) -> Path:
    """Create output directory if needed."""
    path = Path(path_arg).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    return path.resolve()
