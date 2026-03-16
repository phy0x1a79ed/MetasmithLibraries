"""Visualization for functional enrichment results."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)

COMPARISONS = ['POR-S_vs_POR-0', 'POR-8_vs_POR-0', 'POR-8_vs_POR-S']
COMP_LABELS = {
    'POR-S_vs_POR-0': 'Spirulina vs Control',
    'POR-8_vs_POR-0': 'Media #8 vs Control',
    'POR-8_vs_POR-S': 'Media #8 vs Spirulina',
}
NAMESPACE_LABELS = {
    'biological_process': 'Biological Process',
    'molecular_function': 'Molecular Function',
    'cellular_component': 'Cellular Component',
}


def go_dot_plots(go_results, fig_dir, top_n=15):
    """Create GO dot plots: top terms per namespace, sized by count, colored by -log10(padj)."""
    fig_dir = Path(fig_dir)

    for label, df in go_results.items():
        if len(df) == 0:
            continue

        sig = df[df['padj'] < 0.05].copy()
        if len(sig) == 0:
            continue

        sig['neg_log_padj'] = -np.log10(sig['padj'].clip(lower=1e-50))

        namespaces = sig['namespace'].unique()
        n_ns = len(namespaces)
        if n_ns == 0:
            continue

        fig, axes = plt.subplots(1, n_ns, figsize=(7 * n_ns, max(6, min(top_n * 0.4, 10))),
                                  squeeze=False)

        for idx, ns in enumerate(sorted(namespaces)):
            ax = axes[0, idx]
            ns_df = sig[sig['namespace'] == ns].head(top_n).copy()
            if len(ns_df) == 0:
                ax.set_visible(False)
                continue

            ns_df = ns_df.sort_values('neg_log_padj')

            # Truncate long names
            ns_df['display_name'] = ns_df['term_name'].apply(
                lambda x: x[:50] + '...' if len(str(x)) > 50 else x)

            scatter = ax.scatter(
                ns_df['fold_enrichment'],
                range(len(ns_df)),
                s=ns_df['study_count'] * 20,
                c=ns_df['neg_log_padj'],
                cmap='YlOrRd',
                edgecolors='black',
                linewidth=0.5,
                alpha=0.8,
            )
            ax.set_yticks(range(len(ns_df)))
            ax.set_yticklabels(ns_df['display_name'], fontsize=8)
            ax.set_xlabel('Fold Enrichment')
            ax.set_title(NAMESPACE_LABELS.get(ns, ns), fontsize=10)

            plt.colorbar(scatter, ax=ax, label='-log10(padj)', shrink=0.7)

        comp = label.rsplit('_', 1)[0]
        direction = label.rsplit('_', 1)[1]
        fig.suptitle(f'GO Enrichment: {COMP_LABELS.get(comp, comp)} ({direction})',
                     fontsize=12, y=1.02)
        fig.tight_layout()

        for ext in ['pdf', 'png', 'svg']:
            fig.savefig(fig_dir / f'go_dotplot_{label}.{ext}',
                       dpi=300, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"  Saved GO dot plot: {label}")


def kegg_bar_charts(kegg_results, fig_dir, top_n=20):
    """Create KEGG pathway bar charts."""
    fig_dir = Path(fig_dir)

    for label, df in kegg_results.items():
        sig = df[df['padj'] < 0.05].copy() if len(df) > 0 else pd.DataFrame()
        if len(sig) == 0:
            # Show top results even if not significant
            sig = df.head(top_n).copy() if len(df) > 0 else pd.DataFrame()
        if len(sig) == 0:
            continue

        sig = sig.head(top_n).copy()
        sig['neg_log_padj'] = -np.log10(sig['padj'].clip(lower=1e-50))
        sig['display_name'] = sig.apply(
            lambda r: r['pathway_name'][:45] + '...' if len(str(r['pathway_name'])) > 45
            else (r['pathway_name'] if r['pathway_name'] else r['term']),
            axis=1
        )
        sig = sig.sort_values('fold_enrichment')

        fig, ax = plt.subplots(figsize=(8, max(4, len(sig) * 0.35)))

        colors = plt.cm.YlOrRd(sig['neg_log_padj'] / sig['neg_log_padj'].max()
                                if sig['neg_log_padj'].max() > 0 else sig['neg_log_padj'])
        ax.barh(range(len(sig)), sig['fold_enrichment'], color=colors, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(sig)))
        ax.set_yticklabels(sig['display_name'], fontsize=8)
        ax.set_xlabel('Fold Enrichment')

        comp = label.rsplit('_', 1)[0]
        direction = label.rsplit('_', 1)[1]
        ax.set_title(f'KEGG Enrichment: {COMP_LABELS.get(comp, comp)} ({direction})')

        fig.tight_layout()
        for ext in ['pdf', 'png', 'svg']:
            fig.savefig(fig_dir / f'kegg_barplot_{label}.{ext}', dpi=300, bbox_inches='tight')
        plt.close(fig)
        logger.info(f"  Saved KEGG bar chart: {label}")


def cog_grouped_bar(cog_df, fig_dir):
    """Create grouped bar chart of COG category enrichment across comparisons."""
    fig_dir = Path(fig_dir)

    if len(cog_df) == 0:
        return

    # Filter to 'all' direction
    plot_df = cog_df[cog_df['direction'] == 'all'].copy()
    if len(plot_df) == 0:
        return

    # Pivot for grouped bar chart
    pivot = plot_df.pivot_table(
        index='cog_category', columns='comparison',
        values='fold_enrichment', aggfunc='first'
    )

    # Only show categories with data
    pivot = pivot.dropna(how='all')
    if len(pivot) == 0:
        return

    # Add descriptions
    from .cog_analysis import COG_DESCRIPTIONS
    pivot.index = [f"{cat} - {COG_DESCRIPTIONS.get(cat, '')[:30]}" for cat in pivot.index]

    fig, ax = plt.subplots(figsize=(12, max(6, len(pivot) * 0.4)))

    comp_labels = [COMP_LABELS.get(c, c) for c in pivot.columns]
    x = np.arange(len(pivot))
    width = 0.25
    colors = ['#2196F3', '#FF9800', '#4CAF50']

    for i, (col, label) in enumerate(zip(pivot.columns, comp_labels)):
        values = pivot[col].fillna(0)
        ax.barh(x + i * width, values, width, label=label, color=colors[i % 3],
                edgecolor='black', linewidth=0.3)

    ax.axvline(x=1, color='gray', linestyle='--', alpha=0.5, label='Expected')
    ax.set_yticks(x + width)
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xlabel('Fold Enrichment')
    ax.set_title('COG Category Enrichment Across Comparisons')
    ax.legend(fontsize=8)

    fig.tight_layout()
    for ext in ['pdf', 'png', 'svg']:
        fig.savefig(fig_dir / f'cog_grouped_bar.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("  Saved COG grouped bar chart")


def summary_heatmap(go_results, kegg_results, fig_dir, top_n=20):
    """Create summary heatmap of top enriched terms across comparisons."""
    fig_dir = Path(fig_dir)

    # Collect top GO terms across all comparisons (using 'all' direction)
    all_terms = []
    for label, df in go_results.items():
        if not label.endswith('_all'):
            continue
        if len(df) == 0:
            continue
        sig = df[df['padj'] < 0.05].head(top_n)
        for _, row in sig.iterrows():
            all_terms.append({
                'term': row['term'],
                'term_name': row.get('term_name', row['term']),
                'comparison': row['comparison'],
                'neg_log_padj': -np.log10(max(row['padj'], 1e-50)),
                'source': 'GO',
            })

    for label, df in kegg_results.items():
        if not label.endswith('_all'):
            continue
        if len(df) == 0:
            continue
        sig = df[df['padj'] < 0.05].head(10)
        for _, row in sig.iterrows():
            name = row.get('pathway_name', row['term'])
            all_terms.append({
                'term': row['term'],
                'term_name': name if name else row['term'],
                'comparison': row['comparison'],
                'neg_log_padj': -np.log10(max(row['padj'], 1e-50)),
                'source': 'KEGG',
            })

    if not all_terms:
        logger.info("  No significant terms for summary heatmap")
        return

    terms_df = pd.DataFrame(all_terms)

    # Get unique terms (by top significance)
    top_terms = (terms_df.groupby('term_name')['neg_log_padj']
                 .max().sort_values(ascending=False).head(top_n * 2).index)

    # Pivot
    heatmap_data = terms_df[terms_df['term_name'].isin(top_terms)].pivot_table(
        index='term_name', columns='comparison', values='neg_log_padj', aggfunc='max'
    )
    heatmap_data = heatmap_data.fillna(0)

    # Sort by max value
    heatmap_data = heatmap_data.loc[heatmap_data.max(axis=1).sort_values(ascending=False).index]
    heatmap_data = heatmap_data.head(30)

    # Truncate names
    heatmap_data.index = [n[:55] + '...' if len(n) > 55 else n for n in heatmap_data.index]
    heatmap_data.columns = [COMP_LABELS.get(c, c) for c in heatmap_data.columns]

    fig, ax = plt.subplots(figsize=(10, max(6, len(heatmap_data) * 0.35)))
    sns.heatmap(heatmap_data, cmap='YlOrRd', ax=ax, linewidths=0.5,
                cbar_kws={'label': '-log10(padj)'}, annot=True, fmt='.1f')
    ax.set_title('Enrichment Summary: GO + KEGG Terms Across Comparisons')
    ax.set_ylabel('')

    fig.tight_layout()
    for ext in ['pdf', 'png', 'svg']:
        fig.savefig(fig_dir / f'summary_heatmap.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("  Saved summary heatmap")


def mapping_stats_plot(stats, fig_dir):
    """Plot gene ID mapping statistics."""
    fig_dir = Path(fig_dir)

    labels = ['MSTRG\ngenes', 'Mapped to\nBRAKER3', 'With eggNOG\nannotation']
    values = [stats['total_mstrg'], stats['mapped'], stats.get('with_annotation', 0)]
    colors = ['#64B5F6', '#43A047', '#FF9800']

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 50,
                str(val), ha='center', va='bottom', fontweight='bold')

    ax.set_ylabel('Number of Genes')
    ax.set_title('Gene ID Mapping Summary')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    for ext in ['pdf', 'png', 'svg']:
        fig.savefig(fig_dir / f'mapping_stats.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    logger.info("  Saved mapping stats plot")


def create_all_figures(go_results, kegg_results, cog_df, gsea_results, stats, fig_dir):
    """Generate all visualization figures."""
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Creating figures...")

    mapping_stats_plot(stats, fig_dir)
    go_dot_plots(go_results, fig_dir)
    kegg_bar_charts(kegg_results, fig_dir)
    if cog_df is not None and len(cog_df) > 0:
        cog_grouped_bar(cog_df, fig_dir)
    summary_heatmap(go_results, kegg_results, fig_dir)

    logger.info("All figures created.")
