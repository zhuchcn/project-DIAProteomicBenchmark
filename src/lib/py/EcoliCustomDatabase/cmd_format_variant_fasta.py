"""Format variant FASTA file for proteomics analysis."""

from __future__ import annotations
import argparse
from pathlib import Path

from Bio import SeqIO

from Common import setup_logger
from EcoliCustomDatabase.constants import (
    DEFAULT_VARIANT_PEPTIDE_FASTA,
    DEFAULT_FORMATTED_VARIANT_FASTA,
)


LOGGER = setup_logger('EcoliCustomDatabase.FormatVariantFasta')


def build_parser(parser: argparse.ArgumentParser) -> None:
    """Register FormatVariantFasta arguments."""
    parser.description = 'Format variant FASTA headers for proteomics search engines.'
    parser.add_argument(
        '--input-fasta',
        type=Path,
        default=DEFAULT_VARIANT_PEPTIDE_FASTA,
        help=f'Input variant FASTA file (default: {DEFAULT_VARIANT_PEPTIDE_FASTA}).',
    )
    parser.add_argument(
        '--output-fasta',
        type=Path,
        default=DEFAULT_FORMATTED_VARIANT_FASTA,
        help=f'Output formatted FASTA file (default: {DEFAULT_FORMATTED_VARIANT_FASTA}).',
    )


def run(args: argparse.Namespace) -> None:
    """Execute FormatVariantFasta subcommand."""
    input_fasta = _validate_input_file(args.input_fasta, 'Input FASTA')
    output_fasta = _prepare_output_path(args.output_fasta)

    LOGGER.info('Reading input FASTA: %s', input_fasta)
    LOGGER.info('Output FASTA will be written to: %s', output_fasta)

    entries_written = _format_fasta(
        input_fasta,
        output_fasta
    )

    LOGGER.info('Successfully formatted %d entries to: %s', entries_written, output_fasta)


def _format_fasta(input_path: Path, output_path: Path) -> int:
    """
    Format variant FASTA file.

    Split header by spaces to get multiple variant entries

    Returns:
        Number of entries written to output file.
    """
    entries_written = 0

    with output_path.open('w') as outfile:
        for record in SeqIO.parse(input_path, 'fasta'):
            # Original header (without '>')
            header = record.description
            sequence = str(record.seq)

            # Split header by spaces and write separate entries
            entries_written += _write_split_entries(
                outfile,
                header,
                sequence
            )

    return entries_written


def _write_split_entries(outfile, header: str, sequence: str) -> int:
    """
    Split header by spaces and write separate entries for each.

    Returns:
        Number of entries written.
    """
    # Split header by spaces to get individual variant entries
    entries = header.split()

    for entry in entries:
        # Split by '|' to extract fields
        new_header = f'>{entry}'
        outfile.write(new_header + '\n')
        outfile.write(sequence + '\n')

    return len(entries)


def _validate_input_file(path_arg: Path, label: str) -> Path:
    """Validate file existence."""
    path = Path(path_arg).expanduser()
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f'Missing {label}: {path}')
    return path.resolve()


def _prepare_output_path(path_arg: Path) -> Path:
    """Prepare output path and ensure parent directory exists."""
    path = Path(path_arg).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
