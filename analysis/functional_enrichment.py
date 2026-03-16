#!/usr/bin/env python3
"""Functional enrichment analysis for P. purpureum transcriptomics.

Orchestrates gene ID mapping, GO/KEGG ORA, GSEA, COG analysis, and visualization.

Usage:
    python -m analysis.functional_enrichment
    # or
    python analysis/functional_enrichment.py
"""

import sys
import logging
from pathlib import Path

import pandas as pd

# Allow running as script or module
if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from analysis.gene_id_mapping import build_full_mapping, export_clean_mapping
from analysis.enrichment_ora import run_go_enrichment, run_kegg_enrichment
from analysis.enrichment_gsea import run_gsea_analysis
from analysis.cog_analysis import run_cog_analysis
from analysis.visualization import create_all_figures

logger = logging.getLogger(__name__)

REPORT_DIR = 'results/report_package01'
OUTPUT_DIR = 'results/enrichment'
DATA_DIR = 'analysis/data'
FIG_DIR = 'results/enrichment/figures'


def main():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s: %(message)s',
        datefmt='%H:%M:%S',
    )

    report_dir = Path(REPORT_DIR)
    output_dir = Path(OUTPUT_DIR)
    data_dir = Path(DATA_DIR)
    fig_dir = Path(FIG_DIR)

    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Gene ID mapping
    logger.info("=" * 60)
    logger.info("STEP 1: Gene ID Mapping (MSTRG <-> BRAKER3)")
    logger.info("=" * 60)
    mapping_df, eggnog_df, stats = build_full_mapping(report_dir, output_dir)

    # Export clean mapping to report package
    export_clean_mapping(mapping_df, report_dir / 'gene_id_mapping.csv')

    # Load DE results
    de_df = pd.read_csv(report_dir / 'gene_de_results.csv')
    logger.info(f"Loaded {len(de_df)} genes from DE results")

    # Step 2: GO enrichment
    logger.info("=" * 60)
    logger.info("STEP 2: GO Over-Representation Analysis")
    logger.info("=" * 60)
    go_results = run_go_enrichment(de_df, mapping_df, data_dir, output_dir)

    # Step 3: KEGG enrichment
    logger.info("=" * 60)
    logger.info("STEP 3: KEGG Pathway Enrichment")
    logger.info("=" * 60)
    kegg_results = run_kegg_enrichment(de_df, mapping_df, data_dir, output_dir)

    # Step 4: GSEA
    logger.info("=" * 60)
    logger.info("STEP 4: Preranked GSEA")
    logger.info("=" * 60)
    gsea_results = run_gsea_analysis(de_df, mapping_df, data_dir, output_dir)

    # Step 5: COG analysis
    logger.info("=" * 60)
    logger.info("STEP 5: COG Category Analysis")
    logger.info("=" * 60)
    cog_df = run_cog_analysis(de_df, mapping_df, output_dir)

    # Step 6: Visualization
    logger.info("=" * 60)
    logger.info("STEP 6: Visualization")
    logger.info("=" * 60)
    create_all_figures(go_results, kegg_results, cog_df, gsea_results, stats, fig_dir)

    # Summary
    logger.info("=" * 60)
    logger.info("COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Mapping: {stats['mapped']}/{stats['total_mstrg']} MSTRG genes mapped "
                f"({100*stats['mapped']/stats['total_mstrg']:.1f}%)")
    logger.info(f"  With annotations: {stats.get('with_annotation', 0)}")

    n_go_sig = sum(len(df[df['padj'] < 0.05]) for df in go_results.values())
    n_kegg_sig = sum(len(df[df['padj'] < 0.05]) for df in kegg_results.values())
    logger.info(f"GO: {n_go_sig} significant terms across all comparisons")
    logger.info(f"KEGG: {n_kegg_sig} significant pathways across all comparisons")
    logger.info(f"Outputs: {output_dir}/")
    logger.info(f"Figures: {fig_dir}/")


if __name__ == '__main__':
    main()
