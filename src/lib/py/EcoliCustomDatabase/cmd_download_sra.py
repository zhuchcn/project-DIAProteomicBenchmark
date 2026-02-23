"""Download SRA sample and export FASTQ files."""

from __future__ import annotations
import argparse
from pathlib import Path

from Common import DATA_DIR, setup_logger
from EcoliCustomDatabase.constants import DEFAULT_FASTQ_DIR, DEFAULT_SRA_ID, DEFAULT_THREADS
from EcoliCustomDatabase.io_utils import add_conda_env_arg, format_cmd, require_tool, run_cmd


LOGGER = setup_logger('EcoliCustomDatabase.DownloadSRA')
DEFAULT_OUTPUT_DIR = DEFAULT_FASTQ_DIR


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register DownloadSRA arguments."""
    parser.description = 'Download SRA sample and export FASTQ files.'
    add_conda_env_arg(parser)
    parser.add_argument(
        '--sra-id',
        type=str,
        default=DEFAULT_SRA_ID,
        choices=[DEFAULT_SRA_ID],
        help=f'SRA accession to download. Currently fixed to {DEFAULT_SRA_ID}.',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            'Output directory for SRA and FASTQ files. '
            f'Default: {DEFAULT_OUTPUT_DIR}'
        ),
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=DEFAULT_THREADS,
        help=f'Number of threads for fasterq-dump (default: {DEFAULT_THREADS}).',
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Re-run export even if FASTQ output files already exist.',
    )


def run(args: argparse.Namespace) -> None:
    """Execute DownloadSRA subcommand."""
    require_tool('prefetch', args.conda_env)
    require_tool('fasterq-dump', args.conda_env)

    sra_id = args.sra_id
    output_dir = _resolve_output_dir(args.output_dir)
    threads = int(args.threads)
    if threads <= 0:
        raise ValueError(f'--threads must be positive. Got {threads}')

    if _has_fastq_output(output_dir, sra_id) and not args.force:
        LOGGER.info(
            'FASTQ output already exists in %s for %s. Use --force to re-export.',
            output_dir,
            sra_id,
        )
        return

    prefetch_cmd = [
        'prefetch',
        sra_id,
        '--output-directory',
        str(output_dir),
    ]
    LOGGER.info('Running command: %s', format_cmd(prefetch_cmd, args.conda_env))
    run_cmd(prefetch_cmd, args.conda_env)

    sra_path = _find_sra_file(output_dir, sra_id)
    fasterq_cmd = [
        'fasterq-dump',
        str(sra_path),
        '--outdir',
        str(output_dir),
        '--threads',
        str(threads),
    ]
    LOGGER.info('Running command: %s', format_cmd(fasterq_cmd, args.conda_env))
    run_cmd(fasterq_cmd, args.conda_env)

    if not _has_fastq_output(output_dir, sra_id):
        raise RuntimeError(f'FASTQ export appears to have failed for {sra_id} at {output_dir}')

    LOGGER.info('Finished downloading and exporting %s to %s', sra_id, output_dir)


def _resolve_output_dir(path_arg: Path) -> Path:
    """Resolve output directory from DATA_DIR when a relative path is provided."""
    path = path_arg.expanduser()
    if not path.is_absolute():
        path = DATA_DIR / path
    path.mkdir(parents=True, exist_ok=True)
    return path.resolve()


def _find_sra_file(output_dir: Path, sra_id: str) -> Path:
    """Find downloaded SRA file produced by prefetch."""
    candidates = list(output_dir.glob(f'**/{sra_id}.sra'))
    if not candidates:
        raise FileNotFoundError(
            f'Could not find downloaded SRA file for {sra_id} under {output_dir}'
        )
    return candidates[0]


def _has_fastq_output(output_dir: Path, sra_id: str) -> bool:
    """Check whether fasterq-dump output exists."""
    expected = [
        output_dir / f'{sra_id}.fastq',
        output_dir / f'{sra_id}_1.fastq',
        output_dir / f'{sra_id}_2.fastq',
    ]
    return any(path.exists() for path in expected)
