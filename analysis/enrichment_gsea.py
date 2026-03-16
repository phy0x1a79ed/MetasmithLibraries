"""Preranked Gene Set Enrichment Analysis (GSEA)."""

import logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

COMPARISONS = ['POR-S_vs_POR-0', 'POR-8_vs_POR-0', 'POR-8_vs_POR-S']


def build_ranking(de_df, comparison, mapping_df):
    """Build gene ranking: sign(log2FC) * -log10(pvalue).

    Returns: Series indexed by MSTRG gene ID with ranking scores.
    """
    lfc_col = f'{comparison}_log2FoldChange'
    pval_col = f'{comparison}_pvalue'

    # Get genes with valid statistics
    valid = de_df[de_df[pval_col].notna() & de_df[lfc_col].notna()].copy()

    # Only include genes that have annotations
    annotated = set(mapping_df[mapping_df['has_eggnog_annotation']]['mstrg_id'])
    valid = valid[valid['gene_id'].isin(annotated)]

    # Compute ranking metric
    pvals = valid[pval_col].astype(float).clip(lower=1e-300)
    lfcs = valid[lfc_col].astype(float)
    scores = np.sign(lfcs) * -np.log10(pvals)

    ranking = pd.Series(scores.values, index=valid['gene_id'].values)
    ranking = ranking.sort_values(ascending=False)
    # Remove duplicates
    ranking = ranking[~ranking.index.duplicated(keep='first')]

    return ranking


def run_gsea_analysis(de_df, mapping_df, data_dir, output_dir):
    """Run preranked GSEA for GO, KEGG, and COG gene sets."""
    from .enrichment_ora import parse_go_terms, parse_kegg_pathways, load_go_obo, download_go_obo

    output_dir = Path(output_dir)
    data_dir = Path(data_dir)

    # Build gene set collections
    gene_to_go = {}
    gene_to_kegg = {}
    gene_to_cog = {}

    for _, row in mapping_df.iterrows():
        mstrg = row['mstrg_id']

        gos = parse_go_terms(row.get('GOs'))
        if gos:
            gene_to_go[mstrg] = gos

        pathways = parse_kegg_pathways(row.get('KEGG_Pathway'))
        if pathways:
            gene_to_kegg[mstrg] = pathways

        cog = row.get('COG_category')
        if pd.notna(cog) and cog != '-' and cog != '' and cog != 'S':
            gene_to_cog[mstrg] = set(cog)  # Each letter is a category

    # Convert to term -> gene_set format for gseapy
    def invert_mapping(gene_to_terms, min_size=10, max_size=500):
        term_to_genes = defaultdict(set)
        for gene, terms in gene_to_terms.items():
            for term in terms:
                term_to_genes[term].add(gene)
        return {t: list(g) for t, g in term_to_genes.items()
                if min_size <= len(g) <= max_size}

    go_sets = invert_mapping(gene_to_go)
    kegg_sets = invert_mapping(gene_to_kegg, min_size=5, max_size=500)
    cog_sets = invert_mapping(gene_to_cog, min_size=5, max_size=2000)

    logger.info(f"Gene sets — GO: {len(go_sets)}, KEGG: {len(kegg_sets)}, COG: {len(cog_sets)}")

    # Load GO term names for labeling
    obo_path = download_go_obo(data_dir)
    go_term_info = load_go_obo(obo_path)

    all_results = {}

    try:
        import gseapy as gp
    except ImportError:
        logger.warning("gseapy not installed — skipping GSEA analysis")
        return all_results

    for comp in COMPARISONS:
        ranking = build_ranking(de_df, comp, mapping_df)
        logger.info(f"GSEA {comp}: {len(ranking)} ranked genes")

        if len(ranking) < 50:
            logger.info(f"  Skipping — too few ranked genes")
            continue

        for set_name, gene_sets_dict in [('GO', go_sets), ('KEGG', kegg_sets), ('COG', cog_sets)]:
            if not gene_sets_dict:
                continue

            label = f"{comp}_{set_name}"
            logger.info(f"  Running GSEA for {set_name} ({len(gene_sets_dict)} sets)...")

            try:
                result = gp.prerank(
                    rnk=ranking,
                    gene_sets=gene_sets_dict,
                    processes=4,
                    permutation_num=1000,
                    outdir=None,
                    seed=42,
                    min_size=5,
                    max_size=500,
                    verbose=False,
                )

                res_df = result.res2d
                if len(res_df) > 0:
                    res_df = res_df.sort_values('NOM p-val')

                    # Add term names for GO
                    if set_name == 'GO':
                        res_df['term_name'] = res_df['Term'].map(
                            lambda t: go_term_info.get(t, {}).get('name', ''))

                    res_df['comparison'] = comp
                    res_df['gene_set_type'] = set_name

                    out_path = output_dir / f'gsea_{label}.csv'
                    res_df.to_csv(out_path, index=False)

                    n_sig = len(res_df[res_df['FDR q-val'].astype(float) < 0.25])
                    logger.info(f"    {n_sig} significant sets (FDR<0.25)")

                    all_results[label] = res_df

            except Exception as e:
                logger.warning(f"    GSEA failed for {label}: {e}")

    return all_results
