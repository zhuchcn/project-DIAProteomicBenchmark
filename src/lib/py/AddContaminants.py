"""Add cRAP contaminant sequences to a FASTA while cleaning headers."""
from __future__ import annotations
import re
from pathlib import Path

from Template import SubCommand
from Common import (
    setup_logger,
    resolve_file,
)

_PROG_NAME = 'AddContaminants'
_HELP = 'Append cRAP sequences to a FASTA and normalize the headers.'
_DESCRIPTION = (
    'Copy an existing FASTA database, append every entry from the provided '
    'cRAP FASTA, prefix each contaminant header, and strip any PE= annotation.'
)
LOGGER = setup_logger(_PROG_NAME)


class AddContaminants(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--database': {
            'type': Path,
            'required': True,
            'help': 'Path to the existing FASTA database that will be extended.',
        },
        '--crap': {
            'type': Path,
            'required': True,
            'help': 'Path to the cRAP FASTA file whose entries will be appended.',
        },
        '--output': {
            'type': Path,
            'required': True,
            'help': 'Destination path for the expanded FASTA file.',
        },
        '--prefix': {
            'type': str,
            'default': 'contam_',
            'help': 'Prefix to prepend to each contaminant header before the original ID.',
        },
    }

    @staticmethod
    def func(args):
        database_path = resolve_file(args.database, 'FASTA database')
        crap_path = resolve_file(args.crap, 'cRAP FASTA file')
        output_path = args.output.expanduser().resolve()
        prefix = args.prefix or ''

        if not prefix:
            raise ValueError('The contaminant prefix must not be empty.')

        output_path.parent.mkdir(parents=True, exist_ok=True)

        LOGGER.info('Copying base database from %s', database_path)
        last_line_had_newline = True
        with open(database_path, 'r') as base_f, open(output_path, 'w') as out_f:
            for line in base_f:
                out_f.write(line)
                last_line_had_newline = line.endswith('\n')

            if not last_line_had_newline:
                out_f.write('\n')

            LOGGER.info('Appending contaminant entries from %s', crap_path)
            remove_pe = re.compile(r'\s*PE=[^\s]+')
            contaminant_count = 0
            with open(crap_path, 'r') as crap_f:
                for raw_line in crap_f:
                    stripped = raw_line.rstrip('\n')
                    if not stripped:
                        continue
                    if stripped.startswith('>'):
                        cleaned = remove_pe.sub('', stripped[1:]).strip()
                        header_body = f'{prefix}{cleaned}'
                        out_f.write(f'>{header_body}\n')
                        contaminant_count += 1
                    else:
                        out_f.write(f'{stripped}\n')

        LOGGER.info('Added %d contaminant entries; output written to %s', contaminant_count, output_path)
