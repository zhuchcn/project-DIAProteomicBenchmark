"""Annotate variants using VEP."""

from __future__ import annotations
import argparse
from pathlib import Path

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_FILTERED_VCF_GZ,
    DEFAULT_VEP_ASSEMBLY,
    DEFAULT_VEP_CACHE_DIR,
    DEFAULT_VEP_OUTPUT_TSV,
    DEFAULT_VEP_SPECIES,
    DEFAULT_REFERENCE_GENOME
)
from EcoliCustomDatabase.io_utils import add_conda_env_arg, format_cmd, require_tool, run_cmd


LOGGER = setup_logger('EcoliCustomDatabase.AnnotateVEP')


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register AnnotateVEP arguments."""
    parser.description = 'Annotate a VCF using VEP.'
    add_conda_env_arg(parser)
    parser.add_argument(
        '--input-vcf',
        type=Path,
        default=DEFAULT_FILTERED_VCF_GZ,
        help=f'Input VCF/VCF.GZ to annotate (default: {DEFAULT_FILTERED_VCF_GZ}).',
    )
    parser.add_argument(
        '--output-tsv',
        type=Path,
        default=DEFAULT_VEP_OUTPUT_TSV,
        help=f'Output VEP-annotated TSV path (default: {DEFAULT_VEP_OUTPUT_TSV}).',
    )
    parser.add_argument(
        '--species',
        type=str,
        default=DEFAULT_VEP_SPECIES,
        help=f'VEP species (default: {DEFAULT_VEP_SPECIES}).',
    )
    parser.add_argument(
        '--assembly',
        type=str,
        default=DEFAULT_VEP_ASSEMBLY,
        help=f'VEP assembly (default: {DEFAULT_VEP_ASSEMBLY}).',
    )
    parser.add_argument(
        '--dir-cache',
        type=Path,
        default=DEFAULT_VEP_CACHE_DIR,
        help=f'VEP cache directory (default: {DEFAULT_VEP_CACHE_DIR}).',
    )
    parser.add_argument(
        '--fasta',
        type=Path,
        default=DEFAULT_REFERENCE_GENOME,
        help=(
            'Reference genome FASTA file for VEP annotation. '
            f'Default: {DEFAULT_REFERENCE_GENOME}'
        ),
    )
    parser.add_argument(
        '--buffer-size',
        type=int,
        default=10000,
        help='VEP buffer size (default: 10000).',
    )
    parser.add_argument(
        '--distance',
        type=int,
        default=0,
        help='VEP distance value (default: 0).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute AnnotateVEP subcommand."""
    require_tool('vep', args.conda_env)

    input_vcf = _validate_input_file(args.input_vcf, 'input VCF')
    output_tsv = _prepare_output_path(args.output_tsv)
    dir_cache = _prepare_output_dir(args.dir_cache)
    fasta = _validate_input_file(args.fasta, 'reference genome FASTA')

    if args.buffer_size <= 0:
        raise ValueError(f'--buffer-size must be positive. Got {args.buffer_size}')
    if args.distance < 0:
        raise ValueError(f'--distance must be >= 0. Got {args.distance}')

    cmd = [
        'vep',
        '-i',
        str(input_vcf),
        '-o',
        str(output_tsv),
        '--species',
        str(args.species),
        '--assembly',
        str(args.assembly),
        '--cache',
        '--check_ref',
        '--no_stats',
        '--buffer_size',
        str(args.buffer_size),
        '--distance',
        str(args.distance),
        '--dir_cache',
        str(dir_cache),
        '--fasta',
        str(fasta),
        '--offline',
        '--no_intergenic',
        '--force_overwrite',
    ]
    LOGGER.info('Running command: %s', format_cmd(cmd, args.conda_env))
    run_cmd(cmd, args.conda_env)
    LOGGER.info('VEP annotation completed: %s', output_tsv)


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Validate input file existence."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _prepare_output_path(path_arg: Path) -> Path:
    """Prepare output path and create parent directory."""
    path = Path(path_arg).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.resolve()


def _prepare_output_dir(path_arg: Path) -> Path:
    """Prepare output directory."""
    path = Path(path_arg).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    return path.resolve()
