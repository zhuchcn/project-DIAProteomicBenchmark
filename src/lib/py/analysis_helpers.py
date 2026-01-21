""" This module contains functions and class definitions for analysis. """
from __future__ import annotations
from typing import TYPE_CHECKING, Sequence
import re
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pyopenms import (
    PeptideIdentificationList,
    IdXMLFile,
    ProteinIdentification,
    PeptideIdentification,
    PeptideHit,
    AASequence
)


if TYPE_CHECKING:
    from pathlib import Path

TIER_LABEL_MAP = {
    -1: "Contam",
    1: "Ecoli",
    2: "Yeast",
    3: "Human",
}

DEFAULT_TIER_TO_INDIVIDUAL_SAMPLE = {
    "Ecoli": "ecoli_01",
    "Yeast": "yeast_01",
    "Human": "human_01",
}

DEFAULT_TIER_ORDER = ["Ecoli", "Yeast", "Human"]
DEFAULT_STATUS_ORDER = ["individual_only", "mixture_only", "both"]

class IdXMLData:
    def __init__(self, path:Path):
        self.path = path
        self._protein_ids: list[ProteinIdentification] = None
        self._peptide_ids: PeptideIdentificationList = None

    @property
    def protein_ids(self) -> list[ProteinIdentification]:
        if self._protein_ids is None:
            self._load_idxml()
        return self._protein_ids

    @property
    def peptide_ids(self) -> PeptideIdentificationList:
        if self._peptide_ids is None:
            self._load_idxml()
        return self._peptide_ids

    def _load_idxml(self):
        protein_ids = []
        peptide_ids = PeptideIdentificationList()
        IdXMLFile().load(str(self.path), protein_ids, peptide_ids)
        self._protein_ids = protein_ids
        self._peptide_ids = peptide_ids

    def tidy_peptide_id(self, remove_contam:bool = True, remove_decoy: bool = True
            ) -> pd.DataFrame:
        data = {
            'Spectrum Reference': [],
            'RT': [],
            'Observed M/Z': [],
            'Tier': [],
            'Peptide': [],
            'Modified Peptide': [],
            'Unimod Peptide': [],
            'Charge': [],
            'Qvalue': [],
            'Probability': [],
            'Hyperscore': [],
            'Nextscore': []
        }
        pep_id:PeptideIdentification
        for pep_id in self.peptide_ids:
            hits:list[PeptideHit] = pep_id.getHits()
            if not len(hits) > 0:
                continue
            hit:PeptideHit = hits[0]
            if remove_decoy and hit.getMetaValue('target_decoy') == 'decoy':
                continue
            acc = hit.getPeptideEvidences()[0].getProteinAccession()
            if remove_contam and acc.startswith('contam_'):
                continue
            data['Spectrum Reference'].append(pep_id.getMetaValue('spectrum_reference'))
            data['RT'].append(pep_id.getRT())
            data['Observed M/Z'].append(pep_id.getMZ())
            if hit.metaValueExists('tier'):
                tier = TIER_LABEL_MAP[hit.getMetaValue('tier')]
            else:
                tier = None
            data['Tier'].append(tier)
            seq: AASequence = hit.getSequence()
            data['Peptide'].append(seq.toUnmodifiedString())
            data['Modified Peptide'].append(seq.toString())
            data['Unimod Peptide'].append(seq.toUniModString())
            data['Charge'].append(hit.getCharge())
            data['Qvalue'].append(hit.getScore())
            data['Probability'].append(hit.getMetaValue('PeptideProphet probability_score'))
            data['Hyperscore'].append(hit.getMetaValue('hyperscore'))
            data['Nextscore'].append(hit.getMetaValue('nextscore'))
        return pd.DataFrame(data)

    def tidy_protein_groups(self, acc_2_tier, remove_contam: bool = True,
            remove_decoy: bool = True) -> pd.DataFrame:
        df = {
            'Protein Group': [],
            'Probability': [],
            'Qvalue': [],
            'Tier': []
        }
        protein_id = self.protein_ids[0]
        for x, y in zip(protein_id.getProteinGroups(), protein_id.getIndistinguishableProteins()):
            acc = x.accessions[0].decode('utf-8')
            if (acc.startswith('contam_') or acc.startswith('rev_contam_')):
                if remove_contam:
                    continue
                tier = 'Contam'
            else:
                if acc.startswith('rev_'):
                    if remove_decoy:
                        continue
                    acc = acc[4:]
                tier = acc_2_tier[acc]
            df['Protein Group'].append(acc)
            df['Probability'].append(x.probability)
            df['Qvalue'].append(y.probability)
            df['Tier'].append(tier)
        df = pd.DataFrame(df)
        return(df)

class FragPipeResults:
    def __init__(self, base_dir:Path):
        self.base_dir = base_dir
        self._psm = None
        self._peptide = None
        self._protein = None
        self._diann_report = None
        self._split_fdr = None
        self._group_fdr = None

    @property
    def psm(self):
        if self._psm is None:
            self._psm = pd.read_csv(self.base_dir/'psm.tsv', delimiter='\t')
        return self._psm

    @property
    def peptide(self):
        if self._peptide is None:
            self._peptide = pd.read_csv(self.base_dir/'peptide.tsv', delimiter='\t')
        return self._peptide

    @property
    def protein(self):
        if self._protein is None:
            self._protein = pd.read_csv(self.base_dir/'protein.tsv', delimiter='\t')
        return self._protein

    @property
    def diann_report(self):
        if self._protein is None:
            self._protein = pd.read_csv(self.base_dir/'dia-quant-output/report.tsv', delimiter='\t')
        return self._protein

    @property
    def split_fdr(self):
        if self._split_fdr is None:
            base_name = self.base_dir.name
            path = self.base_dir/f'{base_name}.idXML'
            self._split_fdr = IdXMLData(path)
        return self._split_fdr

    @property
    def global_fdr(self):
        if self._group_fdr is None:
            base_name = self.base_dir.name
            path = self.base_dir/f'{base_name}_global_fdr.idXML'
            self._group_fdr = IdXMLData(path)
        return self._group_fdr

def prepare_acc_2_tier_mapping(database_path:Path) -> dict[str, str]:
    acc_2_tier = {}
    tier_map = {
        'PE=1': 'Ecoli',
        'PE=2': 'Yeast',
        'PE=3': 'Human'
    }
    with open(database_path, 'rt') as fh:
        pattern = re.compile(r'PE=\d+')
        for line in fh:
            if not line.startswith('>'):
                continue
            if line.startswith('contam_'):
                continue
            acc = line.split(' ', 1)[0]
            acc = acc.lstrip('>')
            m = pattern.search(line)
            pe = m.group(0) if m else None
            if not pe:
                continue
            tier = tier_map[pe]
            acc_2_tier[acc] = tier
    return acc_2_tier


def compare_mixture_to_individual_protein_groups(
    mixture_tables: dict[str, pd.DataFrame],
    individual_tables: dict[str, pd.DataFrame],
    tier_to_individual_sample: dict[str, str] | None = None,
    tier_column: str = 'Tier',
    protein_column: str = 'Protein Group'
) -> pd.DataFrame:
    """Compile a per-mixture overview that flags cross-sample observation status.

    Each record reports whether the protein group for a given tier is observed in the mixture,
    in the matching individual sample, or in both.
    """
    if tier_to_individual_sample is None:
        tier_to_individual_sample = DEFAULT_TIER_TO_INDIVIDUAL_SAMPLE

    columns = [
        'Mixture Sample',
        'Tier',
        'Protein Group',
        'Individual Sample',
        'Mixture Present',
        'Individual Present',
        'Status'
    ]

    records: list[dict[str, object]] = []
    for mixture_name in sorted(mixture_tables.keys()):
        mixture_df = mixture_tables[mixture_name]
        if mixture_df is None or mixture_df.empty:
            continue
        mixture_df = mixture_df.dropna(subset=[protein_column, tier_column])
        tiers = sorted(mixture_df[tier_column].unique())
        for tier in tiers:
            mixture_groups = set(
                mixture_df.loc[mixture_df[tier_column] == tier, protein_column]
            )
            individual_sample = tier_to_individual_sample.get(tier)
            individual_groups = set()
            if individual_sample:
                individual_df = individual_tables.get(individual_sample)
                if individual_df is not None and not individual_df.empty:
                    individual_df = individual_df.dropna(subset=[protein_column, tier_column])
                    individual_groups = set(
                        individual_df.loc[individual_df[tier_column] == tier, protein_column]
                    )

            all_groups = sorted(mixture_groups | individual_groups)
            for group in all_groups:
                mixture_present = group in mixture_groups
                individual_present = group in individual_groups
                if mixture_present and individual_present:
                    status = 'both'
                elif mixture_present:
                    status = 'mixture_only'
                else:
                    status = 'individual_only'
                records.append({
                    'Mixture Sample': mixture_name,
                    'Tier': tier,
                    'Protein Group': group,
                    'Individual Sample': individual_sample,
                    'Mixture Present': mixture_present,
                    'Individual Present': individual_present,
                    'Status': status
                })

    if not records:
        return pd.DataFrame(columns=columns)
    table = pd.DataFrame(records)
    return table[columns]


def compare_mixture_to_individual_peptides(
    mixture_tables: dict[str, pd.DataFrame],
    individual_tables: dict[str, pd.DataFrame],
    tier_to_individual_sample: dict[str, str] | None = None,
    tier_column: str = 'Tier',
    modified_peptide_column: str = 'Modified Peptide',
    unimod_peptide_column: str = 'Unimod Peptide',
    charge_column: str = 'Charge',
    peptide_column: str = 'Peptide'
) -> pd.DataFrame:
    """Track whether modified peptides appear in mixtures, individuals, or both."""
    if tier_to_individual_sample is None:
        tier_to_individual_sample = DEFAULT_TIER_TO_INDIVIDUAL_SAMPLE

    columns = [
        'Mixture Sample',
        'Tier',
        'Peptide',
        'Modified Peptide',
        'Unimod Peptide',
        'Charge',
        'Individual Sample',
        'Mixture Present',
        'Individual Present',
        'Status'
    ]

    records: list[dict[str, object]] = []
    for mixture_name in sorted(mixture_tables.keys()):
        mixture_df = mixture_tables[mixture_name]
        if mixture_df is None or mixture_df.empty:
            continue
        mixture_df = mixture_df.dropna(subset=[tier_column, modified_peptide_column, charge_column])
        tiers = sorted(mixture_df[tier_column].unique())
        for tier in tiers:
            mixture_subset = mixture_df[mixture_df[tier_column] == tier]
            mixture_keys = set(
                zip(
                    mixture_subset[modified_peptide_column],
                    mixture_subset[charge_column]
                )
            )
            mixture_peptide_map = {
                (row[modified_peptide_column], row[charge_column]):
                    (row.get(peptide_column), row.get(unimod_peptide_column), row.get(charge_column))
                for _, row in mixture_subset.iterrows()
            }

            individual_sample = tier_to_individual_sample.get(tier)
            individual_keys: set[tuple[object, object]] = set()
            individual_peptide_map: dict[tuple[object, object], tuple[object, object, object]] = {}
            if individual_sample:
                individual_df = individual_tables.get(individual_sample)
                if individual_df is not None and not individual_df.empty:
                    individual_df = individual_df.dropna(subset=[modified_peptide_column, charge_column])
                    if tier_column in individual_df.columns and individual_df[tier_column].notna().any():
                        individual_subset = individual_df[individual_df[tier_column] == tier]
                    else:
                        individual_subset = individual_df
                    individual_keys = set(
                        zip(
                            individual_subset[modified_peptide_column],
                            individual_subset[charge_column]
                        )
                    )
                    individual_peptide_map = {
                        (row[modified_peptide_column], row[charge_column]):
                            (row.get(peptide_column), row.get(unimod_peptide_column), row.get(charge_column))
                        for _, row in individual_subset.iterrows()
                    }

            all_keys = sorted(mixture_keys | individual_keys)
            for key in all_keys:
                mixture_present = key in mixture_keys
                individual_present = key in individual_keys
                if mixture_present and individual_present:
                    status = 'both'
                elif mixture_present:
                    status = 'mixture_only'
                else:
                    status = 'individual_only'
                unmod_peptide, unimod_peptide, charge = (
                    mixture_peptide_map.get(key)
                    or individual_peptide_map.get(key)
                )
                records.append({
                    'Mixture Sample': mixture_name,
                    'Tier': tier,
                    'Peptide': unmod_peptide,
                    'Modified Peptide': key[0],
                    'Unimod Peptide': unimod_peptide,
                    'Charge': charge,
                    'Individual Sample': individual_sample,
                    'Mixture Present': mixture_present,
                    'Individual Present': individual_present,
                    'Status': status
                })

    if not records:
        return pd.DataFrame(columns=columns)
    table = pd.DataFrame(records)
    return table[columns]


def summarize_mixture_comparison(
    comparison_table: pd.DataFrame,
    mixture_samples: list[str] | None = None,
    status_column: str = 'Status',
    tier_column: str = 'Tier'
) -> pd.DataFrame:
    """Count proteins per category for each mixture sample and tier."""
    if comparison_table.empty:
        return pd.DataFrame(
            columns=['Mixture Sample', tier_column, status_column, 'Protein Count']
        )
    if mixture_samples is None:
        mixture_samples = sorted(comparison_table['Mixture Sample'].unique())
    filtered = comparison_table[comparison_table['Mixture Sample'].isin(mixture_samples)]
    group_columns = ['Mixture Sample', tier_column, status_column]
    result = (
        filtered
        .groupby(group_columns, dropna=False)
        .size()
        .reset_index(name='Count')
    )
    return result


def plot_status_counts_barplot(
    comparison_table: pd.DataFrame,
    mixture_samples: Sequence[str],
    tier_order: Sequence[str] | None = None,
    status_order: Sequence[str] | None = None,
    count_column: str = 'Count',
    ylabel: str | None = None
) -> plt.Figure:
    """Plot the stacked status counts for each tier and mixture sample."""
    summary = comparison_table
    if count_column not in summary.columns:
        raise ValueError(f'Count column "{count_column}" not found in comparison table.')

    if tier_order is None:
        tier_order = DEFAULT_TIER_ORDER
    if status_order is None:
        status_order = DEFAULT_STATUS_ORDER

    pivot = (
        summary
        .pivot_table(
            index=['Mixture Sample', 'Tier'],
            columns='Status',
            values=count_column,
            fill_value=0
        )
        .reindex(
            pd.MultiIndex.from_product(
                [list(mixture_samples), list(tier_order)],
                names=['Mixture Sample', 'Tier']
            ),
            fill_value=0
        )
    )
    pivot = pivot.reindex(columns=status_order, fill_value=0)

    positions: list[int] = []
    tier_labels: list[str] = []
    counts: dict[str, list[int]] = {status: [] for status in status_order}
    idx = 0
    for sample in mixture_samples:
        for tier in tier_order:
            row = pivot.loc[(sample, tier)]
            tier_labels.append(tier)
            positions.append(idx)
            for status in status_order:
                counts[status].append(row[status])
            idx += 1

    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = [0] * len(positions)
    colors = sns.color_palette('Set2', len(status_order))
    for color, status in zip(colors, status_order):
        ax.bar(
            positions,
            counts[status],
            bottom=bottom,
            label=status.replace('_', ' ').title(),
            color=color
        )
        bottom = [b + v for b, v in zip(bottom, counts[status])]

    ax.set_xticks(positions)
    ax.set_xticklabels(tier_labels)
    ax.set_xlabel('Database Tier')
    ax.set_ylabel(ylabel or count_column.replace('_', ' '))
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    sample_centers: list[float] = []
    tiers_per_sample = len(tier_order)
    for i in range(len(mixture_samples)):
        start = i * tiers_per_sample
        sample_centers.append(start + (tiers_per_sample - 1) / 2)
    ax2.set_xticks(sample_centers)
    ax2.set_xticklabels(mixture_samples)
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.set_xlabel('Mixture Sample')
    ax2.tick_params(axis='x', which='both', length=0)
    ax.legend(title='Status', bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    return fig


def _tier_palette(df: pd.DataFrame) -> dict[str, tuple[float, float, float]]:
    tiers = sorted(df['Tier'].dropna().unique())
    return dict(zip(tiers, sns.color_palette('Set2', n_colors=len(tiers))))


def _log_transform(df: pd.DataFrame, column: str) -> pd.Series:
    return np.log1p(df[column])


def plot_tier_status_grid(df_with_diann: pd.DataFrame) -> plt.Figure:
    """Render the shared grid of Tier-colored plots used in the notebook."""
    panel_df = df_with_diann.copy()
    palette = _tier_palette(panel_df)
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=(1, 2), height_ratios=(2, 1), hspace=0.4, wspace=0.25)

    ax_top_left = fig.add_subplot(gs[0, 0])
    ax_top_right = fig.add_subplot(gs[0, 1], sharey=ax_top_left)
    ax_bottom_left = fig.add_subplot(gs[1, 0])
    ax_bottom_right = fig.add_subplot(gs[1, 1], sharex=ax_top_right)

    def _draw(ax: plt.Axes, subset: pd.DataFrame, column: str, orient: str = 'v') -> None:
        transformed = subset.assign(**{column: _log_transform(subset, column)})
        if orient == 'v':
            sns.boxplot(
                data=transformed,
                x='Tier',
                y=column,
                ax=ax,
                palette=palette,
                showfliers=False,
                whis=[5, 95]
            )
            for tier_name, color in palette.items():
                tier_rows = transformed[transformed['Tier'] == tier_name]
                sns.stripplot(
                    data=tier_rows,
                    x='Tier',
                    y=column,
                    ax=ax,
                    size=4,
                    color=color,
                    jitter=0.2,
                    dodge=True,
                    alpha=0.7
                )
            ax.set_xlabel('Tier')
            ax.set_ylabel(column.replace('_', ' '))
        else:
            sns.boxplot(
                data=transformed,
                x=column,
                y='Tier',
                ax=ax,
                palette=palette,
                orient='h',
                showfliers=False,
                whis=[5, 95]
            )
            for tier_name, color in palette.items():
                tier_rows = transformed[transformed['Tier'] == tier_name]
                sns.stripplot(
                    data=tier_rows,
                    x=column,
                    y='Tier',
                    ax=ax,
                    size=4,
                    color=color,
                    jitter=0.2,
                    dodge=True,
                    alpha=0.7
                )
            ax.set_ylabel('Tier')
            ax.set_xlabel(column.replace('_', ' '))
        ax.grid(True, axis='y' if orient == 'v' else 'x', linestyle='--', alpha=0.6)
        title = column.replace('_', ' ') + f" ({subset['Status'].iloc[0]})" if not subset.empty else 'No data'
        ax.set_title(title)
        if ax.get_legend():
            ax.get_legend().remove()

    status_df = panel_df.dropna(subset=['Tier'])

    mixture_df = status_df[status_df['Status'] == 'mixture_only'].dropna(subset=['PG.MaxLFQ_x'])
    if not mixture_df.empty:
        _draw(ax_top_left, mixture_df, 'PG.MaxLFQ_x', orient='v')
    else:
        ax_top_left.text(0.5, 0.5, 'No mixture_only rows', ha='center', va='center')
        ax_top_left.axis('off')

    both_df = status_df[status_df['Status'] == 'both'].dropna(subset=['PG.MaxLFQ_x', 'PG.MaxLFQ_y'])
    if not both_df.empty:
        both_transformed = both_df.assign(
            PG_MaxLFQ_x=_log_transform(both_df, 'PG.MaxLFQ_x'),
            PG_MaxLFQ_y=_log_transform(both_df, 'PG.MaxLFQ_y')
        )
        sns.scatterplot(
            data=both_transformed,
            x='PG_MaxLFQ_x',
            y='PG_MaxLFQ_y',
            hue='Tier',
            palette=palette,
            ax=ax_top_right,
            edgecolor='black',
            linewidth=0.5
        )
        ax_top_right.set_xlabel('log(PG.MaxLFQ_x + 1)')
        ax_top_right.set_ylabel('log(PG.MaxLFQ_y + 1)')
        ax_top_right.set_title('PG.MaxLFQ (both status)')
        ax_top_right.grid(True, linestyle='--', alpha=0.6)
    else:
        ax_top_right.text(0.5, 0.5, 'No both rows', ha='center', va='center')
        ax_top_right.axis('off')

    ax_bottom_left.axis('off')

    individual_df = status_df[status_df['Status'] == 'individual_only'].dropna(subset=['PG.MaxLFQ_y'])
    if not individual_df.empty:
        _draw(ax_bottom_right, individual_df, 'PG.MaxLFQ_y', orient='h')
    else:
        ax_bottom_right.text(0.5, 0.5, 'No individual_only rows', ha='center', va='center')
        ax_bottom_right.axis('off')

    if not both_df.empty:
        handles, labels = ax_top_right.get_legend_handles_labels()
        ax_top_right.legend(handles, labels, title='Tier', bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax_top_right.legend([], [], frameon=False)

    plt.tight_layout()
    return fig


def plot_status_boxplots(df_with_diann: pd.DataFrame, column_prefix:str) -> plt.Figure:
    """Render a standalone 2×2 grid of Tier-colored PG.MaxLFQ boxplots."""
    panel_df = df_with_diann.copy()
    palette = _tier_palette(panel_df)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    def _draw(ax: plt.Axes, status: str, column: str) -> None:
        subset = panel_df[panel_df['Status'] == status].dropna(subset=[column, 'Tier'])
        if subset.empty:
            ax.text(0.5, 0.5, f'No {status} rows', ha='center', va='center')
            ax.axis('off')
            return
        transformed = subset.assign(**{column: _log_transform(subset, column)})
        sns.boxplot(
            data=transformed,
            x='Tier',
            y=column,
            ax=ax,
            palette=palette,
            showfliers=False,
            whis=[5, 95]
        )
        for tier_name, color in palette.items():
            tier_rows = transformed[transformed['Tier'] == tier_name]
            sns.stripplot(
                data=tier_rows,
                x='Tier',
                y=column,
                ax=ax,
                size=4,
                color=color,
                jitter=0.2,
                dodge=True,
                alpha=0.7
            )
        ax.set_xlabel('Tier')
        ax.set_ylabel(f"log({column} + 1)")
        ax.set_title(f"{column.replace('_', ' ')} ({status})")
        ax.grid(True, axis='y', linestyle='--', alpha=0.6)

    _draw(axes[0, 0], 'mixture_only', f'{column_prefix}_x')
    _draw(axes[0, 1], 'both', f'{column_prefix}_x')
    _draw(axes[1, 0], 'individual_only', f'{column_prefix}_y')
    _draw(axes[1, 1], 'both', f'{column_prefix}_y')

    plt.tight_layout()
    return fig

