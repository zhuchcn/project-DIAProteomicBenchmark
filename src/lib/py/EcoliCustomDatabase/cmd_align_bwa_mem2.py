"""Align paired-end FASTQ reads to bwa-mem2 index."""

from __future__ import annotations
import argparse
from pathlib import Path
import subprocess as sp

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_BWA_MEM2_INDEX_PREFIX,
    DEFAULT_FASTQ_R1,
    DEFAULT_FASTQ_R2,
    DEFAULT_SORTED_BAM,
    DEFAULT_THREADS,
)
from EcoliCustomDatabase.io_utils import (
    add_conda_env_arg,
    format_cmd,
    popen_cmd,
    require_tool,
    run_cmd,
)


LOGGER = setup_logger('EcoliCustomDatabase.AlignBwaMem2')
DEFAULT_READ1 = DEFAULT_FASTQ_R1
DEFAULT_READ2 = DEFAULT_FASTQ_R2
DEFAULT_INDEX_PREFIX = DEFAULT_BWA_MEM2_INDEX_PREFIX
DEFAULT_OUTPUT_BAM = DEFAULT_SORTED_BAM


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register AlignBwaMem2 arguments."""
    parser.description = 'Align reads to reference genome with bwa-mem2.'
    add_conda_env_arg(parser)
    parser.add_argument(
        '--read1',
        type=Path,
        default=DEFAULT_READ1,
        help=f'Path to FASTQ read 1 (default: {DEFAULT_READ1}).',
    )
    parser.add_argument(
        '--read2',
        type=Path,
        default=DEFAULT_READ2,
        help=f'Path to FASTQ read 2 (default: {DEFAULT_READ2}).',
    )
    parser.add_argument(
        '--index-prefix',
        type=Path,
        default=DEFAULT_INDEX_PREFIX,
        help=f'Path prefix of bwa-mem2 index files (default: {DEFAULT_INDEX_PREFIX}).',
    )
    parser.add_argument(
        '--output-bam',
        type=Path,
        default=DEFAULT_OUTPUT_BAM,
        help=f'Output sorted BAM path (default: {DEFAULT_OUTPUT_BAM}).',
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=DEFAULT_THREADS,
        help=f'Number of threads for alignment/sorting (default: {DEFAULT_THREADS}).',
    )
    parser.add_argument(
        '--index-bam',
        action='store_true',
        help='If set, run samtools index on the output BAM.',
    )


def run(args: argparse.Namespace) -> None:
    """Execute AlignBwaMem2 subcommand."""
    require_tool('bwa-mem2', args.conda_env)
    require_tool('samtools', args.conda_env)

    read1 = _validate_input_file(args.read1, 'FASTQ read1')
    read2 = _validate_input_file(args.read2, 'FASTQ read2')
    index_prefix = Path(args.index_prefix).expanduser().resolve()
    output_bam = Path(args.output_bam).expanduser().resolve()
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    threads = int(args.threads)
    if threads <= 0:
        raise ValueError(f'--threads must be positive. Got {threads}')

    _validate_index_prefix(index_prefix)

    mem_cmd = [
        'bwa-mem2',
        'mem',
        '-t',
        str(threads),
        str(index_prefix),
        str(read1),
        str(read2),
    ]
    unsorted_bam = output_bam.with_suffix('.unsorted.bam')
    view_cmd = [
        'samtools',
        'view',
        '-@',
        str(threads),
        '-b',
        '-o',
        str(unsorted_bam),
        '-',
    ]
    sort_cmd = [
        'samtools',
        'sort',
        '-@',
        str(threads),
        '-o',
        str(output_bam),
        str(unsorted_bam),
    ]

    LOGGER.info(
        'Running command: %s | %s',
        format_cmd(mem_cmd, args.conda_env),
        format_cmd(view_cmd, args.conda_env),
    )
    with popen_cmd(mem_cmd, args.conda_env, stdout=sp.PIPE) as mem_proc:
        try:
            run_cmd(view_cmd, args.conda_env, stdin=mem_proc.stdout)
        finally:
            if mem_proc.stdout is not None:
                mem_proc.stdout.close()
        mem_rc = mem_proc.wait()
        if mem_rc != 0:
            raise sp.CalledProcessError(mem_rc, mem_cmd)

    LOGGER.info('Running command: %s', format_cmd(sort_cmd, args.conda_env))
    run_cmd(sort_cmd, args.conda_env)

    if unsorted_bam.exists():
        unsorted_bam.unlink()

    if not output_bam.exists():
        raise RuntimeError(f'Alignment finished but output BAM was not created: {output_bam}')

    if args.index_bam:
        index_cmd = ['samtools', 'index', str(output_bam)]
        LOGGER.info('Running command: %s', format_cmd(index_cmd, args.conda_env))
        run_cmd(index_cmd, args.conda_env)

    LOGGER.info('Alignment completed. Output BAM: %s', output_bam)


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Validate file existence using provided path."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _validate_index_prefix(prefix: Path) -> None:
    """Validate that bwa-mem2 index files exist for the prefix."""
    matches = list(prefix.parent.glob(f'{prefix.name}.*'))
    if not matches:
        raise FileNotFoundError(
            f'No bwa-mem2 index files found for prefix: {prefix}. '
            'Run BuildBwaMem2Index first.'
        )
