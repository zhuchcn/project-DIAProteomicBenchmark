"""Build a training table of FragPipe PSMs annotated with FASTA context."""
from __future__ import annotations
import csv
import re
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from pyopenms import (
	PepXMLFile,
	PeptideIdentificationList,
	PeptideHit,
	PeptideIdentification,
	PeptideEvidence
)

from Template import SubCommand
from Common import resolve_directory, resolve_file, setup_logger

_PROG_NAME = 'PreparePeptideTrainingData'
_HELP = 'Build a peptide-to-protein training table from FragPipe outputs.'
_DESCRIPTION = (
	'Read every pepXML and *_edited.pin file in a FragPipe output directory, ' \
	'merge the scores with FASTA context, and serialize the result as a TSV.'
)
LOGGER = setup_logger(_PROG_NAME)

_PIN_COLUMNS = [
	'unweighted_spectral_entropy',
	'weighted_spectral_entropy',
	'hypergeometric_probability',
	'intersection',
	'pred_RT_real_units',
	'delta_RT_loess'
]

_OUTPUT_COLUMNS = [
	'Peptide',
	'UniModPeptide',
	'Charge',
	'RetentionTime',
	'Probability',
	'HyperScore',
	'ProteinAccession',
	'IsUniquePeptide',
	*(_PIN_COLUMNS),
	'N_Flanking',
	'C_Flanking',
	'N_Term_Cleaved',
	'C_Term_Cleaved',
	'Miscleavages'
]

_PIN_PEPTIDE_PATTERN = re.compile(r'.\.(.+)([0-9])\..$')
_EXPECT_RGX = re.compile(r'(?<=[KR])')

class PreparePeptideTrainingData(SubCommand):
	PROG_NAME = _PROG_NAME
	HELP = _HELP
	DESCRIPTION = _DESCRIPTION
	ARGS = {
		'--fragpipe-dir': {
			'type': Path,
			'required': True,
			'help': 'FragPipe output directory containing pepXML and pin files.'
		},
		'--database': {
			'type': Path,
			'required': True,
			'help': 'FASTA database used in the search.'
		},
		'--output': {
			'type': Path,
			'required': True,
			'help': 'Destination TSV with the merged training table.'
		},
		'--window-size': {
			'type': int,
			'default': 16,
			'help': 'Number of residues to keep for each flanking sequence (default: 16).'
		}
	}

	@staticmethod
	def func(args):
		fragpipe_dir = resolve_directory(args.fragpipe_dir, 'FragPipe output directory')
		pep_files = _collect_pepxml_paths(fragpipe_dir)
		if not pep_files:
			raise FileNotFoundError('No pepXML files were found in the FragPipe directory.')

		pin_files = sorted(fragpipe_dir.glob('*_edited.pin'))
		if not pin_files:
			LOGGER.warning('No *_edited.pin files were detected; spectral entropy fields will be empty.')

		fasta_path = resolve_file(args.database, 'FASTA database file')
		database = _load_fasta(fasta_path)
		metadata = _build_cleavage_map(database)
		canonical_accessions = _build_canonical_accessions(database)
		pin_metrics = _load_pin_metrics(pin_files)
		rows = _collect_training_rows(
			pep_files,
			pin_metrics,
			database,
			metadata,
			canonical_accessions,
			args.window_size
		)

		output_path = args.output.expanduser().resolve()
		output_path.parent.mkdir(parents=True, exist_ok=True)
		with open(output_path, 'w', newline='') as fh:
			writer = csv.DictWriter(fh, fieldnames=_OUTPUT_COLUMNS, delimiter='\t', extrasaction='ignore')
			writer.writeheader()
			writer.writerows(rows)
		LOGGER.info('Wrote %d rows to %s', len(rows), output_path)


def _collect_pepxml_paths(fragpipe_dir: Path) -> list[Path]:
	paths = {*fragpipe_dir.glob('*.pep.xml')}
	return sorted(paths)


def _load_fasta(fasta_path: Path) -> dict[str, str]:
	database = {}
	for record in SeqIO.parse(fasta_path, 'fasta'):
		database[record.id] = str(record.seq)
	return database



def _build_canonical_accessions(database: dict[str, str]) -> dict[str, str]:
	sequence_groups: dict[str, list[str]] = {}
	for accession, sequence in database.items():
		sequence_groups.setdefault(sequence, []).append(accession)
	canonical: dict[str, str] = {}
	for group in sequence_groups.values():
		primary = min(group)
		for accession in group:
			canonical[accession] = primary
	return canonical


def _build_cleavage_map(database: dict[str, str]) -> dict[str, list[int]]:
	cleavage_map = {}
	for acc, seq in database.items():
		sites = [0]
		for match in _EXPECT_RGX.finditer(seq):
			sites.append(match.end())
		sites.append(len(seq))
		cleavage_map[acc] = sorted(set(sites))
	return cleavage_map


def _load_pin_metrics(pin_paths: Iterable[Path]) -> dict[tuple[str, str], dict[str, float | None]]:
	metrics: dict[tuple[str, str], dict[str, float | None]] = {}
	for path in pin_paths:
		with open(path, 'r', newline='') as fh:
			reader = csv.DictReader(fh, delimiter='\t')
			for row in reader:
				cleaned = _clean_pin_peptide(row.get('Peptide', ''))
				if not cleaned:
					continue
				sequence, charge = cleaned
				key = (sequence, charge)
				if key in metrics:
					LOGGER.debug('Skipping duplicate pin entry for %s (charge=%s) in %s', sequence, charge, path.name)
					continue
				metrics[key] = {
					col: _to_float(row.get(col)) for col in _PIN_COLUMNS
				}
	return metrics


def _clean_pin_peptide(peptide: str) -> tuple[str, int] | None:
	match = _PIN_PEPTIDE_PATTERN.search(peptide.strip())
	if not match:
		LOGGER.debug('Unable to parse peptide %s from pin file', peptide)
		return None
	cleaned = match.group(1)
	charge = match.group(2)
	cleaned = cleaned.replace('M[15.9949]', 'M(UniMod:35)')
	cleaned = cleaned.replace('C[57.0215]', 'C(UniMod:4)')
	cleaned = cleaned.replace('n[42.0106]', '.(UniMod:1)')
	return cleaned, int(charge)


def _to_float(value: str | None) -> float | None:
	if value is None:
		return None
	try:
		return float(value)
	except ValueError:
		return None


def _collect_training_rows(
	pep_files: list[Path],
	pin_metrics: dict[tuple[str, str], dict[str, float | None]],
	database: dict[str, str],
	cleavage_map: dict[str, list[int]],
	canonical_accessions: dict[str, str],
	window_size: int
) -> list[dict[str, float | str | bool]]:
	loader = PepXMLFile()
	rows: list[dict[str, float | str | bool]] = []
	for pep_file in pep_files:
		protein_ids = []
		peptide_ids = PeptideIdentificationList()
		loader.load(str(pep_file), protein_ids, peptide_ids)
		identifications = [pep_id for pep_id in peptide_ids]
		LOGGER.info('Processing %d peptide IDs from %s', len(identifications), pep_file.name)
		for pep_id in identifications:
			if not pep_id.getHits():
				continue
			hit = pep_id.getHits()[0]
			if any(_is_decoy_accession(evi.getProteinAccession()) for evi in hit.getPeptideEvidences()):
				continue
			base = _build_base_record(pep_id, hit)
			if base is None:
				continue
			extra_key = (base['UniModPeptide'], base['Charge'])
			pin_data = pin_metrics.get(extra_key)
			target_evidences = [
				evi for evi in hit.getPeptideEvidences()
				if not _is_decoy_accession(evi.getProteinAccession())
			]
			if not target_evidences:
				continue
			canonical_evidence_map: dict[str, PeptideEvidence] = {}
			for evidence in target_evidences:
				accession = evidence.getProteinAccession()
				canonical = canonical_accessions.get(accession, accession)
				canonical_evidence_map.setdefault(canonical, evidence)
			if not canonical_evidence_map:
				continue
			canonical_accessions_sorted = sorted(canonical_evidence_map)
			annotated_rows: list[dict[str, float | str | bool]] = []
			for canonical in canonical_accessions_sorted:
				row = dict(base)
				row['ProteinAccession'] = canonical
				if not _annotate_with_fasta(row, database, cleavage_map, window_size):
					LOGGER.warning(
						'Failed to annotate peptide %s in protein %s; skipping.',
						base['Peptide'], canonical
					)
					continue
				row.update({col: pin_data.get(col) if pin_data else None for col in _PIN_COLUMNS})
				annotated_rows.append(row)
			if not annotated_rows:
				continue
			is_unique = len(annotated_rows) == 1
			for annotated_row in annotated_rows:
				annotated_row['IsUniquePeptide'] = is_unique
				rows.append(annotated_row)
	return rows


def _build_base_record(pep_id: PeptideIdentification, hit: PeptideHit) -> dict[str, float | str]:
	sequence = hit.getSequence()
	if sequence is None:
		return None
	return {
		'Peptide': sequence.toUnmodifiedString(),
		'UniModPeptide': sequence.toUniModString(),
		'Charge': hit.getCharge(),
		'RetentionTime': pep_id.getRT(),
		'Probability': hit.getScore(),
		'HyperScore': _extract_meta_value(hit, ['hyperscore'])
	}


def _extract_meta_value(hit: PeptideHit, keys: list[str]) -> float | str | None:
	for key in keys:
		if hit.metaValueExists(key):
			return hit.getMetaValue(key)
	return None


def _is_decoy_accession(accession: str) -> bool:
	lower = accession.lower()
	return lower.startswith('rev_') or lower.startswith('contam_')


def _annotate_with_fasta(
	row: dict[str, float | str | bool],
	database: dict[str, str],
	cleavage_map: dict[str, list[int]],
	window_size: int
) -> bool:
	protein = row['ProteinAccession']
	sequence = database.get(protein)
	if sequence is None:
		LOGGER.warning('Protein %s not found in the database; skipping row.', protein)
		return False
	peptide = row['Peptide']
	start = sequence.find(peptide)
	if start == -1:
		LOGGER.debug('Peptide %s not found in protein %s; skipping row.', peptide, protein)
		return False
	end = start + len(peptide)
	sites = cleavage_map.get(protein, [])
	row['N_Flanking'] = sequence[max(0, start - window_size):start]
	row['C_Flanking'] = sequence[end:min(len(sequence), end + window_size)]
	row['N_Term_Cleaved'] = 1 if start in sites else 0
	row['C_Term_Cleaved'] = 1 if end in sites else 0
	row['Miscleavages'] = sum(1 for site in sites if start < site < end)
	return True
