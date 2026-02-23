"""Filter proteome FASTA to keep only canonical proteins from Ensembl."""
from __future__ import annotations
from pathlib import Path
from contextlib import ExitStack
from Template import SubCommand
from Common import (
    setup_logger,
    resolve_file,
)
from Bio import SeqIO
from genome_util import parse_gtf


_PROG_NAME = 'FilterCanonicalProteins'
_HELP = 'Filter FASTA to keep only canonical proteins'
_DESCRIPTION = (
    'Filter a proteome FASTA file to keep only proteins from transcripts '
    'marked with the Ensembl_canonical tag in a GTF annotation file.'
)
LOGGER = setup_logger(_PROG_NAME)


def _get_canonical_protein_ids(gtf_path: Path) -> set[str]:
    """
    Parse GTF file and extract canonical protein IDs.

    Args:
        gtf_path: Path to the GTF file

    Returns:
        Set of protein IDs for canonical transcripts
    """
    canonical_protein_ids = set()

    for entry in parse_gtf(gtf_path):
        if 'protein_id' not in entry.attribute:
            continue
        # Check if this transcript has the Ensembl_canonical tag
        tags = entry.attribute.get('tag', [])
        if 'Ensembl_canonical' in tags:
            protein_id = entry.attribute['protein_id']
            canonical_protein_ids.add(protein_id)

    LOGGER.info(f'Found {len(canonical_protein_ids)} canonical proteins')
    return canonical_protein_ids


def _filter_fasta(fasta_path: Path, canonical_protein_ids: set[str], output_path: Path) -> None:
    """
    Filter FASTA file to keep only canonical proteins.

    Args:
        fasta_path: Path to the proteome FASTA file
        canonical_protein_ids: Set of protein IDs to keep
        output_path: Path to write filtered FASTA file
    """
    kept_count = 0
    total_count = 0

    with ExitStack() as stack:
        f_in = stack.enter_context(open(fasta_path, 'r'))
        f_out = stack.enter_context(open(output_path, 'w'))
        writer = SeqIO.FastaIO.FastaWriter(f_out)
        # Parse FASTA file and filter sequences
        for record in SeqIO.parse(f_in, 'fasta'):
            # Extract protein ID (first field separated by space)
            protein_id = record.id.split()[0].split('.')[0]
            total_count += 1

            if protein_id in canonical_protein_ids:
                writer.write_record(record)
                kept_count += 1

    LOGGER.info(f'Processed {total_count} proteins, kept {kept_count} canonical proteins')


class FilterCanonicalProteins(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--gtf': {
            'type': Path,
            'required': True,
            'help': 'Path to the GTF annotation file.',
        },
        '--fasta': {
            'type': Path,
            'required': True,
            'help': 'Path to the proteome FASTA file.',
        },
        '--output': {
            'type': Path,
            'required': True,
            'help': 'Path to the output filtered FASTA file.',
        },
    }

    @staticmethod
    def func(args):
        gtf_path = resolve_file(args.gtf, 'GTF annotation file')
        fasta_path = resolve_file(args.fasta, 'proteome FASTA file')
        output_path = args.output.expanduser().resolve()

        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        LOGGER.info(f'Reading GTF from: {gtf_path}')
        canonical_protein_ids = _get_canonical_protein_ids(gtf_path)

        LOGGER.info(f'Filtering FASTA from: {fasta_path}')
        _filter_fasta(fasta_path, canonical_protein_ids, output_path)

        LOGGER.info(f'Filtered FASTA written to: {output_path}')
