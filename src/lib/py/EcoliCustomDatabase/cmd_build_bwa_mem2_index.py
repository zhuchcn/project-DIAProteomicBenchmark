"""Build bwa-mem2 index for the E. coli reference genome."""

from __future__ import annotations
import argparse
from pathlib import Path

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_BWA_MEM2_INDEX_DIR,
    DEFAULT_BWA_MEM2_INDEX_PREFIX,
    DEFAULT_REFERENCE_GENOME,
)
from EcoliCustomDatabase.io_utils import add_conda_env_arg, format_cmd, require_tool, run_cmd


LOGGER = setup_logger('EcoliCustomDatabase.BuildBwaMem2Index')
DEFAULT_GENOME = DEFAULT_REFERENCE_GENOME
DEFAULT_INDEX_DIR = DEFAULT_BWA_MEM2_INDEX_DIR


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register BuildBwaMem2Index arguments."""
    parser.description = 'Build bwa-mem2 index for the E. coli reference genome.'
    add_conda_env_arg(parser)
    parser.add_argument(
        '--genome',
        type=Path,
        default=DEFAULT_GENOME,
        help=(
            'Reference genome FASTA file. Relative paths are resolved from DATA_DIR. '
            f'Default: {DEFAULT_GENOME}'
        ),
    )
    parser.add_argument(
        '--index-dir',
        type=Path,
        default=DEFAULT_INDEX_DIR,
        help=(
            'Directory to write bwa-mem2 index files. Relative paths are resolved from DATA_DIR. '
            f'Default: {DEFAULT_INDEX_DIR}'
        ),
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default=None,
        help=(
            'Index prefix filename (without path). '
            f'Default: {DEFAULT_BWA_MEM2_INDEX_PREFIX.name}.'
        ),
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='If set, remove existing index files for the target prefix before rebuilding.',
    )


def run(args: argparse.Namespace) -> None:
    """Execute BuildBwaMem2Index subcommand."""
    require_tool('bwa-mem2', args.conda_env)

    genome = _validate_input_file(args.genome, label='reference genome FASTA')
    index_dir = _ensure_output_dir(args.index_dir)
    prefix_name = args.prefix or DEFAULT_BWA_MEM2_INDEX_PREFIX.name
    index_prefix = index_dir / prefix_name

    existing_files = sorted(index_prefix.parent.glob(f'{index_prefix.name}.*'))
    if existing_files and not args.force:
        LOGGER.info(
            'Index files already exist for prefix %s (%d files). Use --force to rebuild.',
            index_prefix,
            len(existing_files),
        )
        return

    if existing_files and args.force:
        for file_path in existing_files:
            file_path.unlink()
        LOGGER.info('Removed %d existing index files for prefix %s', len(existing_files), index_prefix)

    cmd = [
        'bwa-mem2',
        'index',
        '-p',
        str(index_prefix),
        str(genome),
    ]
    LOGGER.info('Running command: %s', format_cmd(cmd, args.conda_env))
    run_cmd(cmd, args.conda_env)

    built_files = sorted(index_prefix.parent.glob(f'{index_prefix.name}.*'))
    if not built_files:
        raise RuntimeError(f'No index files generated for prefix {index_prefix}')
    LOGGER.info('Built bwa-mem2 index at prefix %s (%d files)', index_prefix, len(built_files))


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Ensure provided input file exists without remapping relative paths."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _ensure_output_dir(path_arg: Path) -> Path:
    """Create output directory exactly as provided by caller."""
    path = Path(path_arg).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    return path.resolve()
