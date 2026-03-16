# Functional Enrichment Analysis — *Porphyridium purpureum* Transcriptomics

Differential expression + functional enrichment for *P. purpureum* grown in 3 media conditions (3 replicates each):

| Condition | Description |
|-----------|-------------|
| **POR-0** | Standard *Porphyridium* media (control) |
| **POR-S** | *Spirulina* media formulation |
| **POR-8** | Optimized media #8 — low NaNO₃, low NaCl, high trace metals (Co, Mn, Mo, Zn), high VitB12, high bicarbonate |

## Running

```bash
pip install -r analysis/requirements.txt
python -m analysis.functional_enrichment
```

All outputs go to `results/enrichment/`; figures to `results/enrichment/figures/`.

## Pipeline Summary

### 1. Gene ID Mapping

MSTRG (StringTie) gene IDs mapped to BRAKER3 gene IDs via genomic coordinate overlap using `intervaltree`, then joined to eggNOG annotations.

| Metric | Count |
|--------|-------|
| MSTRG genes (StringTie) | 10,008 |
| BRAKER3 genes | 8,182 |
| Mapped (coordinate overlap) | 8,991 (89.8%) |
| With eggNOG annotation | 6,972 |

**89.8% mapping rate** exceeds the 70% threshold. The 1,017 unmapped MSTRG genes are likely novel transcripts without BRAKER3 gene models (e.g., ncRNAs, UTR-only transcripts).

### 2. GO Over-Representation Analysis (ORA)

Fisher's exact test with BH correction. Background: 3,259 genes with GO annotations.

**171 significant GO terms** across all comparisons (padj < 0.05).

Key downregulated pathways (both POR-S and POR-8 vs control):
- **Ribosome / translation** — cytosolic ribosome (3.0-3.5x), structural constituent of ribosome, cytoplasmic translation. Strong signal in both treatments.
- Organonitrogen compound biosynthesis, amide biosynthesis

Key upregulated pathways:
- **Peroxisome / microbody** (4.3-4.5x enrichment) — upregulated in both POR-S and POR-8 vs control
- Biological regulation, cell communication (POR-8)

### 3. KEGG Pathway Enrichment

Fisher's exact test with BH correction. Background: 2,946 genes with KEGG annotations.

**44 significant KEGG pathways** across all comparisons.

| Comparison | Top Pathway | Fold | padj |
|------------|------------|------|------|
| POR-8 vs POR-0 | Ribosome (ko03010) | 1.9x | 2.3e-04 |
| POR-8 vs POR-0 (down) | Ribosome | 2.6x | 3.8e-08 |
| POR-8 vs POR-0 (down) | Proteasome | 2.9x | 6.1e-03 |
| POR-8 vs POR-0 | Carbon metabolism (ko01200) | 1.7x | 1.1e-02 |
| POR-8 vs POR-S | Alcoholism (histone-related) | 8.8x | 6.5e-05 |

### 4. GSEA (Preranked)

Ranking: `sign(log2FC) × -log10(pvalue)`. 1,000 permutations.

- POR-S vs POR-0: 214 GO sets, 34 KEGG, 5 COG categories (FDR < 0.25)
- POR-8 vs POR-0: 192 GO sets, 29 KEGG, 4 COG categories
- POR-8 vs POR-S: 0 GO, 0 KEGG, 5 COG (fewer DE genes = less power for GO/KEGG, but COG captures broad category shifts)

### 5. COG Category Enrichment

Background: 6,642 genes with COG annotations.

Significant categories (padj < 0.05):

| Category | Description | Enrichment pattern |
|----------|-------------|-------------------|
| **E** | Amino acid transport/metabolism | 1.5x enriched in POR-S DE genes |
| **C** | Energy production/conversion | 1.3-1.4x enriched in both treatments |
| **J** | Translation/ribosomal structure | 1.4x enriched in POR-S |
| **O** | Post-translational modification | 1.2x in POR-S |
| **P** | Inorganic ion transport (trace metals!) | Enriched in DE sets |
| **L** | Replication/recombination/repair | 0.5-0.6x **depleted** in both |
| **T** | Signal transduction | 0.7x depleted in POR-S |
| **B** | Chromatin structure | 3.2x enriched in POR-8 vs POR-S |

## Biological Interpretation

### Consistent across both treatments (vs control):
1. **Translation downregulation** — ribosomal genes massively downregulated. Both media shifts cause translational reprogramming, possibly a stress response reducing growth.
2. **Peroxisome upregulation** — fatty acid β-oxidation, ROS detoxification. Consistent with metabolic stress.
3. **Energy metabolism reshuffling** — COG category C enriched; carbon metabolism KEGG pathway significant.

### POR-8 specific:
- **Proteasome downregulation** — reduced protein turnover
- **Carbon metabolism** enrichment aligns with the high bicarbonate in media #8
- **Chromatin remodeling** (3.2x COG-B enrichment vs POR-S) — epigenetic response to the distinct nutrient profile

### POR-S specific:
- Stronger **amino acid metabolism** signal (COG-E, 1.5x) — may reflect nitrogen source differences between *Spirulina* and *Porphyridium* media

## Challenges & Notes

1. **Gene ID mapping** required coordinate overlap because StringTie (MSTRG.*) and BRAKER3 (g*) use independent gene models. No shared ID system exists between them. Reciprocal overlap with 0.5 threshold works well (89.8% mapping).

2. **Obsolete GO terms** — eggNOG v2.1.12 assigns some obsolete GO terms (e.g., GO:0044445 "cytosolic part"). These appear as unnamed terms in results. They don't affect enrichment validity but cause `NaN` names in a few rows.

3. **POR-8 vs POR-S has few DE genes** (253 with annotations) — ORA finds no significant GO terms. GSEA with COG categories still detects shifts (chromatin, translation) because COG categories are broad enough to capture signal from small gene sets.

4. **KEGG "Alcoholism" and "Systemic lupus"** pathways in POR-8 vs POR-S are false biological leads — these are vertebrate disease pathways that share histone gene families. The real signal is histone/chromatin modification, consistent with the COG-B enrichment.

5. **GO OBO download** requires HTTPS with proper headers (the purl.obolibrary.org redirect returns 403 to bare urllib). Fixed by using a release-specific URL.

## Output Files

```
results/enrichment/
├── gene_id_mapping.csv          # MSTRG <-> BRAKER3 mapping with annotations
├── go_enrichment_*.csv          # GO ORA per comparison × direction
├── kegg_enrichment_*.csv        # KEGG ORA per comparison × direction
├── gsea_*.csv                   # Preranked GSEA results
├── cog_enrichment.csv           # COG Fisher's exact test results
├── cog_summary.csv              # COG × comparison summary table
└── figures/
    ├── mapping_stats.{pdf,png}
    ├── go_dotplot_*.{pdf,png}
    ├── kegg_barplot_*.{pdf,png}
    ├── cog_grouped_bar.{pdf,png}
    └── summary_heatmap.{pdf,png}
```

## Dependencies

```
pandas, numpy, scipy, matplotlib, seaborn, statsmodels
intervaltree    # coordinate overlap mapping
goatools        # GO OBO parsing
gseapy          # preranked GSEA
```
