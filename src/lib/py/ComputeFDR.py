"""Reimagined ComputeFDR workflow (kept separate from the legacy CLI)."""
from __future__ import annotations
from pathlib import Path
import re

from pyopenms import (
    FalseDiscoveryRate,
    IdXMLFile,
    PepXMLFile,
    ProtXMLFile,
    PeptideHit,
    PeptideEvidence,
    PeptideIdentification,
    ProteinHit,
    ProteinIdentification,
    PeptideIdentificationList,
    ProteinGroup,
    String
)

from Template import SubCommand
from Common import (
    resolve_directory,
    resolve_file,
    setup_logger
)

_PROG_NAME = 'ComputeFDRRefactor'
_HELP = 'Refactored ComputeFDR pipeline (new helper-based UX).'
_DESCRIPTION = (
    'Load every pepXML file in a directory, compute peptide/protein FDRs, '
    'and store the resulting identifications to an idXML file.'
)
LOGGER = setup_logger(_PROG_NAME)


class ComputeFDRRefactor(SubCommand):
    PROG_NAME = _PROG_NAME
    HELP = _HELP
    DESCRIPTION = _DESCRIPTION
    ARGS = {
        '--pep': {
            'type': Path,
            'required': True,
            'help': 'Directory containing the pepXML files to process.',
        },
        '--prot': {
            'type': Path,
            'required': False,
            'help': 'Path to the protXML file accompanying the pepXML files.',
        },
        '--output': {
            'type': Path,
            'required': True,
            'help': 'Destination path for the idXML file that will receive the FDR-annotated results.',
        },
        '--decoy-prefix': {
            'type': str,
            'default': 'rev_',
            'help': 'Decoy prefix tag used in the search (default: rev_).'
        },
        '--contam-prefix': {
            'type': str,
            'default': 'contam_',
            'help': 'Contaminant prefix tag used in the search (default: contam_).'
        },
        '--database': {
            'type': Path,
            'default': None,
            'required': False,
            'help': 'FASTA database file used in the search. Required if --group-fdr is set.',
        },
        '--group-fdr': {
            'action': 'store_true',
            'help': 'If set, compute FDRs separately for each database tier indicated '
            "by the PE= annotation in the FASTA headers."
        },
        '--nopg-fdr': {
            'action': 'store_true',
            'help': 'If set, skip protein-group-level FDR computation.'
        }
    }

    @staticmethod
    def func(args):
        pep_dir = resolve_directory(args.pep, 'pepXML directory')
        pep_files = _collect_pepxml_files(pep_dir)
        LOGGER.info('Found %d pepXML files in %s:', len(pep_files), pep_dir)
        if not pep_files:
            raise FileNotFoundError(f'No pepXML files were found in {pep_dir}')

        output_path = args.output.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

        protein_to_tier = _prepare_protein_to_tier_map(args)
        protein_ids, peptide_ids = _load_pepxml_identifications(pep_files)
        if not peptide_ids and not protein_ids:
            raise RuntimeError('No identifications were loaded from the pepXML files.')

        protein_groups, pg_lookup = _resolve_protein_groups(args, protein_ids)

        peptide_ids_tiers = _bucket_peptide_identifications(peptide_ids, args, protein_to_tier)
        protein_ids_tiers = _bucket_protein_hits(protein_ids, args, protein_to_tier)
        protein_group_tiers = _bucket_protein_groups(protein_groups, args, protein_to_tier)

        tiers = _collect_tiers(peptide_ids_tiers, protein_ids_tiers, protein_group_tiers)
        protein_ids_final, peptide_ids_final = _apply_fdr_to_tiers(
            tiers,
            peptide_ids_tiers,
            protein_ids_tiers,
            protein_group_tiers,
            args,
            pg_lookup,
            protein_ids
        )

        IdXMLFile().store(str(output_path), protein_ids_final, peptide_ids_final)
        LOGGER.info('Stored FDR-annotated identifications to %s', output_path)


# ---------------------------------------------------------------------------
# Helper routines
# ---------------------------------------------------------------------------

def _prepare_protein_to_tier_map(args) -> dict[str, int]:
    """Return the protein-to-tier map required for grouped FDR processing."""
    if not args.group_fdr:
        return {}
    if args.database is None:
        raise ValueError('The --database argument must be provided when --group-fdr is set.')
    LOGGER.info('Generating protein to tier mapping from database %s', args.database)
    database_fasta = resolve_file(args.database, 'FASTA database file')
    return _generate_protein_to_tier_map(database_fasta)


def _load_pepxml_identifications(
    pep_files: list[Path]
) -> tuple[list[ProteinIdentification], PeptideIdentificationList]:
    """Load protein and peptide identifications from every pepXML file."""
    loader = PepXMLFile()
    protein_ids: list[ProteinIdentification] = []
    peptide_ids = PeptideIdentificationList()
    for pep_file in pep_files:
        LOGGER.info('Loading pepXML %s', pep_file)
        prot_ids = []
        pep_ids = PeptideIdentificationList()
        loader.load(str(pep_file), prot_ids, pep_ids)
        protein_ids.extend(prot_ids)
        for pep in pep_ids:
            peptide_ids.push_back(pep)
    return protein_ids, peptide_ids


def _resolve_protein_groups(
    args,
    protein_ids: list[ProteinIdentification]
) -> tuple[ProteinIdentification, dict]:
    """Return the resolved protein groups plus a lookup for probability propagation."""
    if args.nopg_fdr:
        if args.prot:
            LOGGER.warning('The --nopg-fdr flag is set, ignoring the provided protXML file.')
        return ProteinIdentification(), {}

    prot_file = resolve_file(args.prot, 'protXML file') if args.prot else None
    protein_groups = ProteinIdentification()
    if prot_file is not None:
        placeholder = PeptideIdentification()
        LOGGER.info('Loading protein groups from protXML %s', prot_file)
        ProtXMLFile().load(str(prot_file), protein_groups, placeholder)
    else:
        for prot_id in protein_ids:
            for pg in prot_id.getProteinGroups():
                protein_groups.insertProteinGroup(pg)
    pg_lookup = {pg.accessions[0]: pg for pg in protein_groups.getProteinGroups()}
    return protein_groups, pg_lookup


def _bucket_peptide_identifications(
    peptide_ids: PeptideIdentificationList,
    args,
    protein_to_tier: dict[str, int]
) -> dict[float, PeptideIdentificationList]:
    """Organize peptide-level identifications into tiered buckets."""
    peptide_tiers: dict[float, PeptideIdentificationList] = {}
    for pep in peptide_ids:
        pep = _infer_psm_db_info(
            peptide_id=pep,
            contam_prefix=args.contam_prefix,
            decoy_prefix=args.decoy_prefix,
            group_fdr=args.group_fdr,
            protein_to_tier=protein_to_tier
        )
        tier = pep.getHits()[0].getMetaValue('tier') if args.group_fdr and pep.getHits() else 0
        peptide_tiers.setdefault(tier, PeptideIdentificationList()).push_back(pep)
    return peptide_tiers


def _bucket_protein_hits(
    protein_ids: list[ProteinIdentification],
    args,
    protein_to_tier: dict[str, int]
) -> dict[float, list[ProteinHit]]:
    """Assign protein hits to tier buckets for independent FDR application."""
    protein_tiers: dict[float, list[ProteinHit]] = {}
    for protein_id in protein_ids:
        for hit in protein_id.getHits():
            hit = _infer_protein_db_info(
                protein_hit=hit,
                contam_prefix=args.contam_prefix,
                decoy_prefix=args.decoy_prefix,
                group_fdr=args.group_fdr,
                protein_to_tier=protein_to_tier
            )
            tier = hit.getMetaValue('tier') if args.group_fdr else 0
            protein_tiers.setdefault(tier, []).append(hit)
    return protein_tiers


def _bucket_protein_groups(
    protein_groups: ProteinIdentification,
    args,
    protein_to_tier: dict[str, int]
) -> dict[float, list[ProteinGroup]]:
    """Organize protein groups into tiered buckets (even when empty)."""
    group_tiers: dict[float, list[ProteinGroup]] = {}
    if args.nopg_fdr:
        return group_tiers
    for pg in protein_groups.getProteinGroups():
        tier, pg_filtered = _infer_protein_group_db_info(
            protein_group=pg,
            contam_prefix=args.contam_prefix,
            decoy_prefix=args.decoy_prefix,
            group_fdr=args.group_fdr,
            protein_to_tier=protein_to_tier
        )
        group_tiers.setdefault(tier, []).append(pg_filtered)
    return group_tiers


def _collect_tiers(
    *tier_maps: dict[float, object]
) -> set[float]:
    """Return the union of keys from all tier dictionaries."""
    tiers: set[float] = set()
    for tier_map in tier_maps:
        tiers.update(tier_map.keys())
    return tiers


def _apply_fdr_to_tiers(
    tiers: set[float],
    peptide_ids_tiers: dict[float, PeptideIdentificationList],
    protein_ids_tiers: dict[float, list[ProteinHit]],
    protein_group_tiers: dict[float, list[ProteinGroup]],
    args,
    pg_lookup: dict,
    protein_ids: list[ProteinIdentification]
) -> tuple[list[ProteinIdentification], PeptideIdentificationList]:
    """Execute peptide- and protein-level FDR per tier and merge the results."""
    peptide_ids_final = PeptideIdentificationList()
    protein_ids_final: list[ProteinIdentification] = []
    protein_id_final = ProteinIdentification()

    fdr = FalseDiscoveryRate()
    params = fdr.getParameters()
    params.setValue("add_decoy_proteins", "true")
    fdr.setParameters(params)

    for tier in sorted(tiers):
        pep_ids = PeptideIdentificationList()
        for pep in peptide_ids_tiers.get(tier, []):
            pep_ids.push_back(pep)

        prot_hits = protein_ids_tiers.get(tier, [])
        prot_ids = ProteinIdentification()
        prot_ids.setHits(prot_hits)
        if tier in protein_group_tiers:
            for pg in protein_group_tiers[tier]:
                prot_ids.insertProteinGroup(pg)
                prot_ids.insertIndistinguishableProteins(pg)

        LOGGER.info('Processing tier %s: %d peptide IDs, %d protein IDs', tier, len(pep_ids), len(prot_hits))
        fdr.apply(pep_ids)
        if args.nopg_fdr:
            fdr.applyBasic(prot_ids, groups_too=False)
        else:
            fdr.applyPickedProteinFDR(
                id=prot_ids,
                decoy_string=String(args.decoy_prefix),
                decoy_prefix=True,
                groups_too=True
            )

        for pep in pep_ids:
            peptide_ids_final.push_back(pep)
        for prot_hit in prot_ids.getHits():
            protein_id_final.insertHit(prot_hit)

        if not args.nopg_fdr:
            for pg_i in prot_ids.getProteinGroups():
                pg_o = pg_lookup[pg_i.accessions[0]]
                pg_o.probability = pg_i.probability
                protein_id_final.insertProteinGroup(pg_o)
            for ip_i in prot_ids.getIndistinguishableProteins():
                ip_o = pg_lookup[ip_i.accessions[0]]
                ip_o.probability = ip_i.probability
                protein_id_final.insertIndistinguishableProteins(ip_o)

    if not protein_ids:
        raise RuntimeError('Unable to set protein metadata; no ProteinIdentification entries are available.')
    reference_protein = protein_ids[0]
    protein_id_final.setIdentifier(reference_protein.getIdentifier())
    protein_id_final.setHigherScoreBetter(reference_protein.isHigherScoreBetter())
    ms_runs = []
    reference_protein.getPrimaryMSRunPath(ms_runs)
    protein_id_final.setPrimaryMSRunPath(ms_runs)
    protein_id_final.setScoreType(reference_protein.getScoreType())
    protein_id_final.setSearchEngine(reference_protein.getSearchEngine())
    protein_id_final.setSearchEngineVersion(reference_protein.getSearchEngineVersion())
    protein_id_final.setSearchParameters(reference_protein.getSearchParameters())
    protein_id_final.setSignificanceThreshold(reference_protein.getSignificanceThreshold())
    protein_id_final.setDateTime(reference_protein.getDateTime())

    protein_ids_final.append(protein_id_final)
    return protein_ids_final, peptide_ids_final


# ---------------------------------------------------------------------------
# Existing inference helpers
# ---------------------------------------------------------------------------

def _infer_psm_db_info(peptide_id: PeptideIdentification, contam_prefix:str,
        decoy_prefix:str, group_fdr:bool, protein_to_tier:dict[str, int]) -> PeptideIdentification:
    """ Infer the database information from the PSMs in the provided PeptideIdentification.

    The peptide sequence of each `PeptideHit` can have multiple `PeptideEvidence`,
    because a peptide can map to multiple proteins. These proteins are sorted such that
    contaminants are on the top, then decoys, then ordered by tier if `group_fdr` is True.
    """
    hits = []
    hit:PeptideHit
    for hit in peptide_id.getHits():
        evi: PeptideEvidence
        evidences = []
        for evi in hit.getPeptideEvidences():
            acc = evi.getProteinAccession()
            is_contam = acc.startswith(contam_prefix) or acc.startswith(decoy_prefix + contam_prefix)
            is_target = not acc.startswith(decoy_prefix)
            if group_fdr:
                if is_contam:
                    tier = -1
                else:
                    tier = protein_to_tier.get(acc, float('inf'))
                evidences.append((evi, is_target, tier))
            else:
                evidences.append((evi, is_target, 0))
        evidences_sorted = sorted(evidences,  key=lambda x: (x[1], x[2]))
        is_target = evidences_sorted[0][1]
        target_decoy = 'target' if is_target else 'decoy'
        hit.setMetaValue('target_decoy', target_decoy)
        if group_fdr:
            tier = evidences_sorted[0][2]
            hit.setMetaValue('tier', tier)
        hit.setPeptideEvidences([evi[0] for evi in evidences_sorted])
        hits.append(hit)
    # Are PeptideHits always sorted that the best score is first?
    peptide_id.setHits(hits)
    return peptide_id


def _infer_protein_db_info(protein_hit: ProteinHit, contam_prefix:str,
        decoy_prefix:str, group_fdr:bool, protein_to_tier:dict[str, int]) -> ProteinHit:
    """ Infer the database information from the provided ProteinHit. """
    acc = protein_hit.getAccession()
    is_contam = acc.startswith(contam_prefix) or acc.startswith(decoy_prefix + contam_prefix)
    is_target = not acc.startswith(decoy_prefix)
    if group_fdr:
        if is_contam:
            tier = -1
        else:
            tier = protein_to_tier.get(acc, float('inf'))
    else:
        tier = 0
    target_decoy = 'target' if is_target else 'decoy'
    protein_hit.setMetaValue('target_decoy', target_decoy)
    protein_hit.setMetaValue('tier', tier)
    return protein_hit


def _infer_protein_group_db_info(protein_group: ProteinGroup, contam_prefix:str,
        decoy_prefix:str, group_fdr:bool, protein_to_tier:dict[str, int]
        ) -> tuple[int, ProteinGroup]:
    """ Infer the database information from the provided ProteinGroup. """
    if not group_fdr:
        return 0, protein_group
    tiers = []
    for acc in protein_group.accessions:
        acc = acc.decode('utf-8')
        is_contam = acc.startswith(contam_prefix) or acc.startswith(decoy_prefix + contam_prefix)
        if is_contam:
            tier = -1
        else:
            tier = protein_to_tier.get(acc, float('inf'))
        tiers.append(tier)
    tier = tiers[0]
    protein_group_filtered = ProteinGroup()
    protein_group_filtered.probability = protein_group.probability
    protein_group_filtered.accessions = [
        acc_i for acc_i, tier_i in zip(protein_group.accessions, tiers) if tier_i == tier
    ]
    return tier, protein_group_filtered


def _collect_pepxml_files(directory: Path) -> list[Path]:
    """Return every pepXML file contained in the provided directory."""
    pattern = '*.pep.xml'
    return sorted(directory.glob(pattern))


def _generate_protein_to_tier_map(fasta_path:Path) -> dict[str, int]:
    """Generate a mapping from peptide sequence to its tier based on the provided FASTA database."""
    protein_to_tier = {}
    pattern = re.compile(r'PE=(\d+)')
    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            if not line.startswith('>'):
                continue
            m = pattern.search(line)
            if not m:
                continue
            tier = int(m.group(1))
            protein = line.split(' ', 1)[0][1:]
            protein_to_tier[protein] = tier
    return protein_to_tier
