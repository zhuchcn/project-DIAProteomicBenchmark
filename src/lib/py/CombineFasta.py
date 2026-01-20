"""Combine multiple FASTA files into one with PE= annotations."""
from __future__ import annotations
import csv
import re
from pathlib import Path

from Template import SubCommand
from Common import (
	setup_logger,
	resolve_file,
)


_PROG_NAME = 'CombineFasta'
_HELP = 'Combine multiple FASTA files with PE= annotations'
_DESCRIPTION = (
	'Read a TSV manifest file and combine multiple FASTA files into one, '
	'adding PE= keyword to each header with the database order value.'
)
LOGGER = setup_logger(_PROG_NAME)


class CombineFasta(SubCommand):
	PROG_NAME = _PROG_NAME
	HELP = _HELP
	DESCRIPTION = _DESCRIPTION
	ARGS = {
		'--manifest': {
			'type': Path,
			'required': True,
			'help': (
				'Path to the TSV manifest file with 3 columns: '
				'database_name, database_order, fasta_path'
			),
		},
		'--output': {
			'type': Path,
			'required': True,
			'help': 'Path to the output combined FASTA file.',
		},
	}

	@staticmethod
	def func(args):
		manifest_path = resolve_file(args.manifest, 'manifest TSV file')
		output_path = args.output.expanduser().resolve()

		# Ensure output directory exists
		output_path.parent.mkdir(parents=True, exist_ok=True)

		LOGGER.info(f'Reading manifest from: {manifest_path}')

		# Read the manifest file
		with open(manifest_path, 'r') as f:
			reader = csv.DictReader(f, delimiter='\t')
			manifest_entries = list(reader)

		if not manifest_entries:
			raise ValueError('Manifest file is empty')

		# Validate that the required columns exist
		required_columns = {'database_name', 'database_order', 'fasta_path'}
		actual_columns = set(manifest_entries[0].keys())

		if not required_columns.issubset(actual_columns):
			raise ValueError(
				f'Manifest must contain columns: {required_columns}. '
				f'Found: {actual_columns}'
			)

		LOGGER.info(f'Found {len(manifest_entries)} FASTA files to combine')

		validated_entries = []
		for entry in manifest_entries:
			db_name = entry['database_name']
			db_order = entry['database_order']
			fasta_path = Path(entry['fasta_path']).expanduser()
			try:
				resolved_path = resolve_file(fasta_path, f'FASTA file for {db_name}')
			except FileNotFoundError as e:
				LOGGER.error(str(e))
				raise
			validated_entries.append((db_name, db_order, resolved_path))

		LOGGER.info(f'Validated {len(validated_entries)} FASTA files exist')

		# Open output file for writing
		with open(output_path, 'w') as out_f:
			for db_name, db_order, fasta_path in validated_entries:
				LOGGER.info(
					f'Processing {db_name} (order={db_order}) from {fasta_path}'
				)

				# Read and process the FASTA file
				with open(fasta_path, 'r') as in_f:
					for line in in_f:
						line = line.rstrip('\n')
						if line.startswith('>'):
							# This is a header line - remove any existing PE= entries
							cleaned_header = re.sub(r'\s*PE=[^\s]+', '', line).rstrip()
							modified_line = f'{cleaned_header} PE={db_order}'
							out_f.write(modified_line + '\n')
						else:
							# This is a sequence line - write as-is
							out_f.write(line + '\n')

		LOGGER.info(f'Successfully created combined FASTA file: {output_path}')
