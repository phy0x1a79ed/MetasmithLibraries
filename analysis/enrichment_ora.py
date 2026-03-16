"""Over-representation analysis (ORA) for GO terms and KEGG pathways."""

import logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

COMPARISONS = ['POR-S_vs_POR-0', 'POR-8_vs_POR-0', 'POR-8_vs_POR-S']


def parse_go_terms(go_str):
    """Parse GO terms from eggNOG GOs column."""
    if pd.isna(go_str) or go_str == '-' or go_str == '':
        return set()
    return set(go_str.split(','))


def parse_kegg_pathways(kegg_str):
    """Parse KEGG pathways, keeping only ko* prefixed entries."""
    if pd.isna(kegg_str) or kegg_str == '-' or kegg_str == '':
        return set()
    pathways = set()
    for p in kegg_str.split(','):
        p = p.strip()
        if p.startswith('ko'):
            pathways.add(p)
    return pathways


def parse_kegg_ko(ko_str):
    """Parse KEGG KO terms."""
    if pd.isna(ko_str) or ko_str == '-' or ko_str == '':
        return set()
    return set(ko_str.replace('ko:', '').split(','))


def get_de_gene_sets(de_df, comparison, padj_cutoff=0.05, lfc_cutoff=1.0):
    """Extract DE gene sets for a comparison.

    Returns dict with keys: 'all', 'up', 'down' -> set of gene_ids
    """
    padj_col = f'{comparison}_padj'
    lfc_col = f'{comparison}_log2FoldChange'

    if padj_col not in de_df.columns:
        return {'all': set(), 'up': set(), 'down': set()}

    sig = de_df[de_df[padj_col].notna() & (de_df[padj_col] < padj_cutoff)].copy()
    sig_lfc = sig[sig[lfc_col].abs() >= lfc_cutoff]

    up = set(sig_lfc[sig_lfc[lfc_col] > 0]['gene_id'])
    down = set(sig_lfc[sig_lfc[lfc_col] < 0]['gene_id'])

    return {'all': up | down, 'up': up, 'down': down}


def run_fisher_enrichment(study_genes, background_genes, gene_to_terms, min_term_size=3):
    """Run Fisher's exact test for term enrichment.

    Args:
        study_genes: set of DE gene IDs
        background_genes: set of all gene IDs with annotations
        gene_to_terms: dict mapping gene_id -> set of terms
        min_term_size: minimum genes per term to test

    Returns: DataFrame with enrichment results
    """
    # Build term -> gene sets
    term_to_genes = defaultdict(set)
    for gene, terms in gene_to_terms.items():
        if gene in background_genes:
            for term in terms:
                term_to_genes[term].add(gene)

    study_in_bg = study_genes & background_genes
    n_bg = len(background_genes)
    n_study = len(study_in_bg)

    results = []
    for term, term_genes in term_to_genes.items():
        if len(term_genes) < min_term_size:
            continue

        # 2x2 contingency table
        a = len(study_in_bg & term_genes)  # DE and in term
        b = len(study_in_bg - term_genes)  # DE and not in term
        c = len(term_genes - study_in_bg)  # not DE and in term
        d = n_bg - a - b - c              # not DE and not in term

        if a == 0:
            continue

        _, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        fold_enrichment = (a / n_study) / (len(term_genes) / n_bg) if n_study > 0 else 0

        results.append({
            'term': term,
            'study_count': a,
            'study_total': n_study,
            'bg_count': len(term_genes),
            'bg_total': n_bg,
            'fold_enrichment': fold_enrichment,
            'pvalue': pval,
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    df = df.sort_values('pvalue')

    # BH correction
    if len(df) > 0:
        _, padj, _, _ = multipletests(df['pvalue'], method='fdr_bh')
        df['padj'] = padj

    return df


def load_go_obo(obo_path):
    """Parse go-basic.obo to get GO term names and namespaces."""
    terms = {}
    current_id = None
    current_name = None
    current_namespace = None
    is_obsolete = False

    with open(obo_path) as f:
        for line in f:
            line = line.strip()
            if line == '[Term]':
                if current_id and not is_obsolete:
                    terms[current_id] = {'name': current_name, 'namespace': current_namespace}
                current_id = None
                current_name = None
                current_namespace = None
                is_obsolete = False
            elif line.startswith('id: GO:'):
                current_id = line.split(': ', 1)[1]
            elif line.startswith('name: '):
                current_name = line.split(': ', 1)[1]
            elif line.startswith('namespace: '):
                current_namespace = line.split(': ', 1)[1]
            elif line == 'is_obsolete: true':
                is_obsolete = True

    # Don't forget last term
    if current_id and not is_obsolete:
        terms[current_id] = {'name': current_name, 'namespace': current_namespace}

    return terms


def download_go_obo(data_dir):
    """Download go-basic.obo if not present."""
    import urllib.request
    obo_path = Path(data_dir) / 'go-basic.obo'
    if not obo_path.exists():
        logger.info("Downloading go-basic.obo...")
        url = 'https://release.geneontology.org/2024-06-17/ontology/go-basic.obo'
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                with open(obo_path, 'wb') as f:
                    f.write(resp.read())
            logger.info(f"  Saved to {obo_path}")
        except Exception:
            # Fallback: try goatools downloader
            logger.info("  Primary download failed, trying goatools...")
            from goatools.base import download_go_basic_obo
            download_go_basic_obo(str(obo_path))
    return obo_path


def fetch_kegg_pathway_names(data_dir):
    """Fetch KEGG pathway names from REST API or cache."""
    import urllib.request
    cache_path = Path(data_dir) / 'kegg_pathway_names.tsv'
    if cache_path.exists():
        names = {}
        with open(cache_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    names[parts[0]] = parts[1]
        return names

    logger.info("Fetching KEGG pathway names...")
    try:
        url = 'https://rest.kegg.jp/list/pathway/ko'
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read().decode('utf-8')

        names = {}
        with open(cache_path, 'w') as f:
            for line in data.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    # path:ko00010 -> ko00010
                    pathway_id = parts[0].replace('path:', '')
                    name = parts[1]
                    names[pathway_id] = name
                    f.write(f"{pathway_id}\t{name}\n")

        logger.info(f"  Cached {len(names)} pathway names")
        return names
    except Exception as e:
        logger.warning(f"Could not fetch KEGG pathway names: {e}")
        return {}


def run_go_enrichment(de_df, mapping_df, data_dir, output_dir):
    """Run GO over-representation analysis for all comparisons."""
    output_dir = Path(output_dir)
    data_dir = Path(data_dir)

    # Download/load GO OBO
    obo_path = download_go_obo(data_dir)
    go_terms = load_go_obo(obo_path)
    logger.info(f"Loaded {len(go_terms)} GO terms from OBO")

    # Build gene-to-GO mapping (using MSTRG IDs)
    gene_to_go = {}
    for _, row in mapping_df.iterrows():
        gos = parse_go_terms(row.get('GOs'))
        if gos:
            gene_to_go[row['mstrg_id']] = gos

    background = set(gene_to_go.keys())
    logger.info(f"GO background: {len(background)} genes with GO annotations")

    all_results = {}
    for comp in COMPARISONS:
        gene_sets = get_de_gene_sets(de_df, comp)
        for direction, genes in gene_sets.items():
            label = f"{comp}_{direction}"
            mapped_genes = genes & background
            logger.info(f"  {label}: {len(mapped_genes)} DE genes with GO annotations "
                        f"(of {len(genes)} total DE)")

            if len(mapped_genes) < 5:
                logger.info(f"    Skipping — too few genes")
                continue

            result = run_fisher_enrichment(mapped_genes, background, gene_to_go)
            if len(result) > 0:
                # Add GO term names and namespaces
                result['term_name'] = result['term'].map(
                    lambda t: go_terms.get(t, {}).get('name', ''))
                result['namespace'] = result['term'].map(
                    lambda t: go_terms.get(t, {}).get('namespace', ''))
                result['comparison'] = comp
                result['direction'] = direction

                out_path = output_dir / f'go_enrichment_{label}.csv'
                result.to_csv(out_path, index=False)
                logger.info(f"    {len(result[result['padj'] < 0.05])} significant terms (padj<0.05)")

                all_results[label] = result

    return all_results


def run_kegg_enrichment(de_df, mapping_df, data_dir, output_dir):
    """Run KEGG pathway over-representation analysis."""
    output_dir = Path(output_dir)
    data_dir = Path(data_dir)

    # Fetch pathway names
    pathway_names = fetch_kegg_pathway_names(data_dir)

    # Build gene-to-KEGG mapping
    gene_to_kegg = {}
    for _, row in mapping_df.iterrows():
        pathways = parse_kegg_pathways(row.get('KEGG_Pathway'))
        if pathways:
            gene_to_kegg[row['mstrg_id']] = pathways

    background = set(gene_to_kegg.keys())
    logger.info(f"KEGG background: {len(background)} genes with KEGG annotations")

    all_results = {}
    for comp in COMPARISONS:
        gene_sets = get_de_gene_sets(de_df, comp)
        for direction, genes in gene_sets.items():
            label = f"{comp}_{direction}"
            mapped_genes = genes & background

            if len(mapped_genes) < 5:
                continue

            result = run_fisher_enrichment(mapped_genes, background, gene_to_kegg, min_term_size=2)
            if len(result) > 0:
                result['pathway_name'] = result['term'].map(
                    lambda t: pathway_names.get(t, ''))
                result['comparison'] = comp
                result['direction'] = direction

                out_path = output_dir / f'kegg_enrichment_{label}.csv'
                result.to_csv(out_path, index=False)
                logger.info(f"  {label}: {len(result[result['padj'] < 0.05])} significant pathways")

                all_results[label] = result

    return all_results
