"""COG category enrichment analysis."""

import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

COMPARISONS = ['POR-S_vs_POR-0', 'POR-8_vs_POR-0', 'POR-8_vs_POR-S']

COG_DESCRIPTIONS = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Post-translational modification, chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'X': 'Mobilome: prophages, transposons',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton',
}


def get_de_gene_sets(de_df, comparison, padj_cutoff=0.05, lfc_cutoff=1.0):
    padj_col = f'{comparison}_padj'
    lfc_col = f'{comparison}_log2FoldChange'

    if padj_col not in de_df.columns:
        return {'all': set(), 'up': set(), 'down': set()}

    sig = de_df[de_df[padj_col].notna() & (de_df[padj_col] < padj_cutoff)].copy()
    sig_lfc = sig[sig[lfc_col].abs() >= lfc_cutoff]
    up = set(sig_lfc[sig_lfc[lfc_col] > 0]['gene_id'])
    down = set(sig_lfc[sig_lfc[lfc_col] < 0]['gene_id'])
    return {'all': up | down, 'up': up, 'down': down}


def run_cog_analysis(de_df, mapping_df, output_dir):
    """Analyze COG category distribution and enrichment."""
    output_dir = Path(output_dir)

    # Build gene-to-COG mapping
    gene_to_cog = {}
    for _, row in mapping_df.iterrows():
        cog = row.get('COG_category')
        if pd.notna(cog) and cog != '-' and cog != '':
            # Multi-letter COG = gene belongs to multiple categories
            gene_to_cog[row['mstrg_id']] = set(cog)

    background = set(gene_to_cog.keys())
    logger.info(f"COG background: {len(background)} genes with COG annotations")

    # Count background distribution
    bg_counts = defaultdict(int)
    for gene, cats in gene_to_cog.items():
        for cat in cats:
            bg_counts[cat] += 1

    all_results = []

    for comp in COMPARISONS:
        gene_sets = get_de_gene_sets(de_df, comp)

        for direction in ['all', 'up', 'down']:
            de_genes = gene_sets[direction] & background
            if len(de_genes) < 5:
                continue

            # Count DE genes per category
            de_counts = defaultdict(int)
            for gene in de_genes:
                for cat in gene_to_cog.get(gene, set()):
                    de_counts[cat] += 1

            n_de = len(de_genes)
            n_bg = len(background)

            for cat in sorted(set(bg_counts.keys()) | set(de_counts.keys())):
                a = de_counts.get(cat, 0)
                b = n_de - a
                c = bg_counts.get(cat, 0) - a
                d = n_bg - a - b - c

                if a == 0:
                    continue

                _, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
                expected = n_de * bg_counts.get(cat, 0) / n_bg
                fold = a / expected if expected > 0 else 0

                all_results.append({
                    'comparison': comp,
                    'direction': direction,
                    'cog_category': cat,
                    'cog_description': COG_DESCRIPTIONS.get(cat, 'Unknown'),
                    'de_count': a,
                    'de_total': n_de,
                    'bg_count': bg_counts.get(cat, 0),
                    'bg_total': n_bg,
                    'fold_enrichment': fold,
                    'pvalue': pval,
                })

    if not all_results:
        return pd.DataFrame()

    df = pd.DataFrame(all_results)

    # BH correction per comparison+direction group
    padj_values = []
    for (comp, direction), group in df.groupby(['comparison', 'direction']):
        if len(group) > 0:
            _, padj, _, _ = multipletests(group['pvalue'], method='fdr_bh')
            padj_values.extend(zip(group.index, padj))

    padj_series = pd.Series(dict(padj_values))
    df['padj'] = padj_series

    out_path = output_dir / 'cog_enrichment.csv'
    df.to_csv(out_path, index=False)
    logger.info(f"Saved COG analysis to {out_path}")

    # Summary table: COG × comparison showing fold enrichment
    summary_rows = []
    for cat in sorted(COG_DESCRIPTIONS.keys()):
        row = {'cog_category': cat, 'description': COG_DESCRIPTIONS[cat]}
        for comp in COMPARISONS:
            mask = (df['comparison'] == comp) & (df['direction'] == 'all') & (df['cog_category'] == cat)
            subset = df[mask]
            if len(subset) > 0:
                row[f'{comp}_fold'] = subset.iloc[0]['fold_enrichment']
                row[f'{comp}_padj'] = subset.iloc[0]['padj']
            else:
                row[f'{comp}_fold'] = np.nan
                row[f'{comp}_padj'] = np.nan
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dir / 'cog_summary.csv', index=False)

    return df
