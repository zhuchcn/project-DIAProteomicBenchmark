"""Call and filter variants from aligned BAM using bcftools."""

from __future__ import annotations
import argparse
from pathlib import Path
import subprocess as sp

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_FILTERED_VCF_GZ,
    DEFAULT_RAW_VCF_GZ,
    DEFAULT_REFERENCE_GENOME,
    DEFAULT_SORTED_BAM,
)
from EcoliCustomDatabase.io_utils import (
    add_conda_env_arg,
    format_cmd,
    popen_cmd,
    require_tool,
    run_cmd,
)


LOGGER = setup_logger('EcoliCustomDatabase.CallVariantBcftools')
DEFAULT_FILTER_EXPR = 'QUAL<30 || DP<10'


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register CallVariantBcftools arguments."""
    parser.description = 'Call and filter variants using bcftools.'
    add_conda_env_arg(parser)
    parser.add_argument(
        '--genome',
        type=Path,
        default=DEFAULT_REFERENCE_GENOME,
        help=f'Reference genome FASTA path (default: {DEFAULT_REFERENCE_GENOME}).',
    )
    parser.add_argument(
        '--bam',
        type=Path,
        default=DEFAULT_SORTED_BAM,
        help=f'Input aligned/sorted BAM path (default: {DEFAULT_SORTED_BAM}).',
    )
    parser.add_argument(
        '--output-vcf',
        type=Path,
        default=DEFAULT_RAW_VCF_GZ,
        help=f'Output raw VCF.GZ path from mpileup+call (default: {DEFAULT_RAW_VCF_GZ}).',
    )
    parser.add_argument(
        '--output-filtered-vcf',
        type=Path,
        default=DEFAULT_FILTERED_VCF_GZ,
        help=(
            'Output filtered VCF.GZ path from bcftools filter '
            f'(default: {DEFAULT_FILTERED_VCF_GZ}).'
        ),
    )
    parser.add_argument(
        '--filter-expr',
        type=str,
        default=DEFAULT_FILTER_EXPR,
        help=f'bcftools filter expression (default: {DEFAULT_FILTER_EXPR!r}).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute CallVariantBcftools subcommand."""
    require_tool('bcftools', args.conda_env)

    genome = _validate_input_file(args.genome, 'reference genome FASTA')
    bam = _validate_input_file(args.bam, 'input BAM')
    output_vcf = _prepare_output_path(args.output_vcf)
    output_filtered_vcf = _prepare_output_path(args.output_filtered_vcf)
    filter_expr = str(args.filter_expr)

    mpileup_cmd = [
        'bcftools',
        'mpileup',
        '-f',
        str(genome),
        str(bam),
    ]
    call_cmd = [
        'bcftools',
        'call',
        '-mv',
        '-Oz',
        '-o',
        str(output_vcf),
    ]

    LOGGER.info(
        'Running command: %s | %s',
        format_cmd(mpileup_cmd, args.conda_env),
        format_cmd(call_cmd, args.conda_env),
    )
    with popen_cmd(mpileup_cmd, args.conda_env, stdout=sp.PIPE) as mpileup_proc:
        try:
            run_cmd(call_cmd, args.conda_env, stdin=mpileup_proc.stdout)
        finally:
            if mpileup_proc.stdout is not None:
                mpileup_proc.stdout.close()
        mpileup_rc = mpileup_proc.wait()
        if mpileup_rc != 0:
            raise sp.CalledProcessError(mpileup_rc, mpileup_cmd)

    if not output_vcf.exists():
        raise RuntimeError(f'Raw VCF output was not created: {output_vcf}')

    filter_cmd = [
        'bcftools',
        'filter',
        '-e',
        filter_expr,
        str(output_vcf),
        '-Oz',
        '-o',
        str(output_filtered_vcf),
    ]
    LOGGER.info('Running command: %s', format_cmd(filter_cmd, args.conda_env))
    run_cmd(filter_cmd, args.conda_env)

    if not output_filtered_vcf.exists():
        raise RuntimeError(f'Filtered VCF output was not created: {output_filtered_vcf}')

    LOGGER.info(
        'Variant calling completed. Raw VCF: %s, Filtered VCF: %s',
        output_vcf,
        output_filtered_vcf,
    )


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Validate file existence using provided path."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _prepare_output_path(path_arg: Path) -> Path:
    """Prepare output path and create parent directory."""
    path = Path(path_arg).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.resolve()
