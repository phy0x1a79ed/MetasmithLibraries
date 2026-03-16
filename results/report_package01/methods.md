# Methods

## Metabolomics Processing and Pathway Analysis

### Workflow orchestration

All bioinformatics steps from raw LC-MS data to pathway enrichment were orchestrated using Metasmith (https://github.com/hallamlab/Metasmith), a workflow planning engine that resolves computational dependencies as a directed acyclic graph (DAG) and dispatches execution via Nextflow (Di Tommaso et al. 2017). Each computational step ran in an isolated container environment managed by Apptainer (formerly Singularity; Kurtzer et al. 2017). The complete pipeline comprised 5 computational steps (Figure 1).

### Data acquisition

Targeted and untargeted LC-MS metabolomics data were acquired for *Porphyridium purpureum* strain 161 (161-PE) under 4 conditions: baseline (0h), selenium stress (8), phosphorus deprivation (P), and sulfur deprivation (S), each measured at 3 timepoints (24h, 72h, 144h) with 4 biological replicates per condition-timepoint. Two chromatographic methods (C18 reverse-phase and HILIC) were used, each in positive and negative ionization modes, yielding 4 targeted datasets and 4 untargeted datasets.

Targeted analysis used JGI's Metatlas workflow with compound atlas databases for compound identification. Untargeted feature detection was performed by MZmine with GNPS2 Feature-Based Molecular Networking (FBMN) for library matching.

### Normalization

Targeted peak heights were normalized using Internal Standard (ISTD) correction. Per-sample normalization factors were computed as the median ISTD peak height divided by the global median of medians. Raw intensities were divided by these factors to correct for injection volume and matrix effects.

Untargeted features were filtered (3x extraction control threshold) and only features with GNPS2 library matches were retained. Targeted and untargeted features were merged, with duplicates (same InChI key and polarity) resolved in favor of targeted identifications.

All intensities were log2-transformed with zero imputation (min positive value / 2).

### Compound annotation

Targeted compounds were annotated from JGI CompoundAtlas databases providing KEGG IDs, InChI keys, PubChem CIDs, HMDB IDs, ChEBI IDs, and LipidMaps IDs. Untargeted compounds were annotated by GNPS2-FBMN spectral library matching. Compounds lacking KEGG IDs but possessing PubChem CIDs were queried against the KEGG REST API (conv/compound/pubchem) for KEGG ID conversion.

### Differential analysis

Differential metabolite abundance was computed using Welch's t-test (scipy.stats.ttest_ind, equal_var=False) for 9 pairwise comparisons (3 conditions x 3 timepoints vs 0h baseline). Log2 fold changes were computed directly from log2-transformed means. P-values were corrected using Benjamini-Hochberg FDR correction. Metabolites were considered significant at adjusted p-value < 0.05 and |log2 fold change| > 1.0.

One-way ANOVA (scipy.stats.f_oneway) was also performed per timepoint to test for differences across all 3 treatment conditions simultaneously.

### Pathway enrichment

Over-Representation Analysis (ORA) was performed using Fisher's exact test (one-sided, greater alternative) for each KEGG pathway in each comparison. The background set comprised all detected metabolites with KEGG compound IDs (n = varies per analysis). Pathway-compound mappings were retrieved from the KEGG REST API (Kanehisa et al. 2023). Pathways with fewer than 2 background compounds were excluded. P-values were corrected using Benjamini-Hochberg FDR.

### Visualization

Volcano plots (3x3 grid, conditions x timepoints) were generated using matplotlib v3.10.8 (Hunter 2007). Pathway enrichment was visualized as dot plots with size proportional to overlap ratio and color mapped to -log10(adjusted p-value). The pipeline DAG was rendered using the graphviz Python package v0.21 (Ellson et al. 2001). All figures were exported in SVG (vector) and PNG (300 DPI) formats.

Clustering: Gaussian Mixture Model (GMM) clustering was applied to a combined 18-feature vector per metabolite (log2FC + -log10(p-value) × 3 timepoints × 3 conditions), standardized with StandardScaler. The optimal number of components was selected via BIC elbow criterion (k=2..9, first k where marginal improvement drops below 5%). Cluster labels (upregulated/downregulated/insignificant) were assigned based on mean log2FC of each component centroid (threshold ±0.3).

Dimensionality reduction for visualization used supervised UMAP (McInnes et al. 2018) with GMM cluster labels as the target variable (n_neighbors=10, min_dist=0.8, spread=2.0, repulsion_strength=1.5, target_metric="categorical"). A 4-panel UMAP figure shows GMM clusters (top-left) and per-nutrient fold-change/significance overlays (remaining panels). A 3×3 cluster-colored volcano grid (conditions × timepoints) was also generated.

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.x | Analysis runtime |
| pandas | 2.3.3 | Data manipulation |
| numpy | 2.4.3 | Numerical operations |
| scipy | 1.16.3 | Statistical tests (Welch's t-test, Fisher's exact, ANOVA) |
| matplotlib | 3.10.8 | Visualization |
| graphviz | 0.21 | DAG rendering |
| requests | — | KEGG REST API access |
| scikit-learn | — | GMM clustering, StandardScaler |
| umap-learn | — | Supervised UMAP dimensionality reduction |

### Data summary

| Metric | Value |
|--------|-------|
| Targeted compounds | 1642 (after dedup) |
| Sample columns | 160 |
| Comparisons | 9 |
| Chromatographies | C18, HILIC |
| Polarities | POS, NEG |

---

## References

Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289-300.

Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316-319. https://doi.org/10.1038/nbt.3820

Ellson, J., Gansner, E., Koutsofios, L., North, S. & Woodhull, G. (2001). Graphviz — open source graph drawing tools. In *International Symposium on Graph Drawing*, pp. 483-484.

Hunter, J.D. (2007). Matplotlib: a 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90-95. https://doi.org/10.1109/MCSE.2007.55

Kanehisa, M. et al. (2023). KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Research*, 51(D1), D587-D592. https://doi.org/10.1093/nar/gkac963

Kurtzer, G.M., Sochat, V. & Bauer, M.W. (2017). Singularity: scientific containers for mobility of compute. *PLOS ONE*, 12(5), e0177459. https://doi.org/10.1371/journal.pone.0177459

Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261-272. https://doi.org/10.1038/s41592-019-0686-2

---

*Analysis performed by Claude (Opus 4.6, Anthropic)*
*2026-03-16*
