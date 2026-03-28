#!/usr/bin/env python3
"""Reconcile organellar vs nuclear transcriptomics results and compile report_package02.

Loads pipeline outputs (organellar counts, DE results, contamination report),
identifies assembly contigs matching organellar sequences, removes contaminated
MSTRG genes from nuclear tables, and builds a reconciled report package.
"""

import csv
import logging
import shutil
from pathlib import Path

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

# Paths
MLIB = Path(__file__).resolve().parent.parent
RESULTS_BASE = Path.home() / "agentic_workspace/results/metasmith-libraries/transcriptomics/pipeline"
REPORT_PKG01 = RESULTS_BASE / "report_package01"
REPORT_PKG02 = RESULTS_BASE / "report_package02"

# Pipeline output locations (adjust after pipeline run if needed)
PIPELINE_RESULTS = RESULTS_BASE / "organellar_pipeline_results"

BLAST_HEADER = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]


def load_contamination_report(path):
    """Load BLASTn contamination report and identify contaminated contigs."""
    if not path.exists() or path.stat().st_size == 0:
        logger.info("No contamination detected (empty report)")
        return pd.DataFrame(columns=BLAST_HEADER), set()

    df = pd.read_csv(path, sep="\t", header=None, names=BLAST_HEADER)
    logger.info(f"Loaded {len(df)} BLAST hits from contamination report")

    # Identify assembly contigs with significant organellar matches
    # Filter: >=90% identity and alignment length >=500bp (likely full organellar contigs)
    significant = df[(df["pident"] >= 90) & (df["length"] >= 500)]
    contaminated_contigs = set(significant["sseqid"].unique())

    # Also flag shorter matches (potential NUMTs/NUPTs) but don't auto-remove
    short_matches = df[(df["pident"] >= 90) & (df["length"] >= 100) & (df["length"] < 500)]
    numt_contigs = set(short_matches["sseqid"].unique()) - contaminated_contigs

    if contaminated_contigs:
        logger.info(f"Contaminated contigs (>=500bp match): {contaminated_contigs}")
    if numt_contigs:
        logger.info(f"Potential NUMT/NUPT contigs (100-500bp match, flagged only): {numt_contigs}")

    return df, contaminated_contigs


def load_gene_id_mapping(path):
    """Load gene_id_mapping.csv from report_package01."""
    df = pd.read_csv(path)
    logger.info(f"Loaded gene ID mapping: {len(df)} genes")
    return df


def identify_organellar_mstrg_genes(mapping_df, contaminated_contigs):
    """Find MSTRG genes located on organellar-contaminated contigs."""
    if not contaminated_contigs:
        return set()

    organellar_genes = set(
        mapping_df[mapping_df["contig"].isin(contaminated_contigs)]["mstrg_id"]
    )
    logger.info(f"Found {len(organellar_genes)} MSTRG genes on contaminated contigs")
    return organellar_genes


def build_organellar_gene_mapping(blast_df, mapping_df, contaminated_contigs):
    """Map MSTRG genes on contaminated contigs to organellar genes via BLAST coordinates."""
    if blast_df.empty or not contaminated_contigs:
        return pd.DataFrame()

    # For each contaminated contig, get BLAST alignment coordinates
    # Then match MSTRG genes on that contig to organellar genes by coordinate overlap
    rows = []
    contam_genes = mapping_df[mapping_df["contig"].isin(contaminated_contigs)]

    for _, gene in contam_genes.iterrows():
        mstrg_id = gene["mstrg_id"]
        contig = gene["contig"]
        m_start = gene["mstrg_start"]
        m_end = gene["mstrg_end"]

        # Find BLAST hits on this contig
        contig_hits = blast_df[blast_df["sseqid"] == contig]
        for _, hit in contig_hits.iterrows():
            # Check if MSTRG gene overlaps with the BLAST alignment on the assembly
            s_start = min(hit["sstart"], hit["send"])
            s_end = max(hit["sstart"], hit["send"])

            overlap_start = max(m_start, s_start)
            overlap_end = min(m_end, s_end)

            if overlap_start <= overlap_end:
                overlap_len = overlap_end - overlap_start + 1
                gene_len = m_end - m_start + 1
                overlap_frac = overlap_len / gene_len if gene_len > 0 else 0

                rows.append({
                    "mstrg_id": mstrg_id,
                    "contig": contig,
                    "mstrg_start": m_start,
                    "mstrg_end": m_end,
                    "organellar_seq": hit["qseqid"],
                    "organellar_start": hit["qstart"],
                    "organellar_end": hit["qend"],
                    "blast_pident": hit["pident"],
                    "blast_length": hit["length"],
                    "overlap_fraction": round(overlap_frac, 4),
                })

    result = pd.DataFrame(rows)
    if not result.empty:
        # Keep best match per MSTRG gene
        result = result.sort_values("overlap_fraction", ascending=False).drop_duplicates(
            subset=["mstrg_id"], keep="first"
        )
    logger.info(f"Built organellar gene mapping: {len(result)} MSTRG→organellar matches")
    return result


def clean_nuclear_tables(count_path, de_path, organellar_mstrg_genes, output_dir):
    """Remove organellar genes from nuclear count and DE tables."""
    # Count table
    counts = pd.read_csv(count_path, index_col="gene_id")
    n_before = len(counts)
    counts_clean = counts[~counts.index.isin(organellar_mstrg_genes)]
    n_removed = n_before - len(counts_clean)
    logger.info(f"Nuclear count table: {n_before} → {len(counts_clean)} genes ({n_removed} removed)")
    counts_clean.to_csv(output_dir / "nuclear_gene_count_table.csv")

    # DE results
    de = pd.read_csv(de_path, index_col="gene_id")
    de_clean = de[~de.index.isin(organellar_mstrg_genes)]
    logger.info(f"Nuclear DE results: {len(de)} → {len(de_clean)} genes")
    de_clean.to_csv(output_dir / "nuclear_de_results.csv")

    return counts_clean, de_clean


def compile_methods(output_dir, n_organellar_genes, n_contam_contigs, n_removed):
    """Write methods description for report_package02."""
    methods = f"""\
# Organellar Transcriptomics Methods

## Organellar Reference Extraction

Organellar genome sequences and gene annotations were extracted from NCBI GenBank
reference records: NC_023133.1 (chloroplast, ~217 kb) and MT483997.1 (mitochondrion,
~129 kb) for *Porphyridium purpureum*. BioPython `SeqIO.parse()` was used to extract
FASTA sequences and convert feature annotations to GFF3 format, preserving gene names,
products, locus tags, and translation tables (transl_table=4 for mitochondrial genes
where TGA encodes tryptophan).

## Prokaryotic RNA-seq Alignment

RNA-seq reads from all 9 samples (POR-0/POR-S/POR-8 × 3 replicates) were aligned to
each organellar genome using minimap2 v2.28 in short-read mode (`-x sr`), which does
not model splice junctions — appropriate for prokaryotic-style gene organization in
organellar genomes. Alignments were sorted and indexed with samtools.

## Organellar Gene Quantification

Read counts per organellar gene were computed using pysam, counting primary alignments
overlapping gene intervals from the GFF3 annotations. A gene-level count matrix was
generated with gene IDs, gene names, and per-sample counts.

## Organellar Differential Expression

PyDESeq2 was used for differential expression analysis of organellar genes, following
the same pairwise comparison design as the nuclear analysis (POR-S vs POR-0, POR-8 vs
POR-0, POR-8 vs POR-S). Significance thresholds: adjusted p-value < 0.05, |log2FC| > 1.

## Contamination Assessment

BLASTn was used to search organellar sequences against the chromosome-level genome
assembly to identify potential organellar contamination. Hits with ≥90% identity and
≥500 bp alignment length were classified as contaminated contigs. Shorter matches
(100-500 bp) were flagged as potential NUMTs/NUPTs but not removed.

## Nuclear-Organellar Reconciliation

{n_contam_contigs} assembly contig(s) were identified as organellar contamination,
containing {n_removed} MSTRG gene(s). These genes were removed from the nuclear
count and DE tables. The organellar pipeline provides independent expression
estimates for {n_organellar_genes} organellar genes aligned directly to the
reference organellar genomes.
"""
    with open(output_dir / "methods.md", "w") as f:
        f.write(methods)
    logger.info("Wrote methods.md")


def compile_results_narrative(output_dir, org_counts, org_de, blast_df,
                              contaminated_contigs, n_removed, org_mapping):
    """Write narrative results for report_package02."""
    # Summarize organellar expression
    n_org_genes = len(org_counts) if org_counts is not None else 0

    # Count DEGs per comparison
    deg_summary = []
    if org_de is not None:
        for col in org_de.columns:
            if col.endswith("_padj"):
                comparison = col.replace("_padj", "")
                lfc_col = f"{comparison}_log2FoldChange"
                if lfc_col in org_de.columns:
                    sig = (org_de[col].fillna(1) < 0.05) & (org_de[lfc_col].abs() > 1)
                    deg_summary.append(f"- {comparison}: {sig.sum()} DEGs")

    narrative = f"""\
# Organellar Transcriptomics Analysis Results

## Overview

Organellar transcriptomics analysis was performed on *Porphyridium purpureum*
chloroplast and mitochondrial genomes, covering {n_org_genes} organellar genes
across 9 samples (POR-0, POR-S, POR-8 conditions, 3 replicates each).

## Contamination Assessment

{"No significant organellar contamination was detected in the chromosome assembly."
 if not contaminated_contigs else
 f"{len(contaminated_contigs)} assembly contig(s) were identified as containing organellar sequences: {', '.join(sorted(contaminated_contigs))}. "
 f"{n_removed} MSTRG gene(s) on these contigs were removed from the nuclear tables."}

{f"Total BLAST hits: {len(blast_df)}" if not blast_df.empty else ""}

## Organellar Differential Expression

{chr(10).join(deg_summary) if deg_summary else "No differential expression results available."}

## Reconciliation Summary

- Nuclear genes before cleanup: reported in report_package01
- Organellar genes removed from nuclear tables: {n_removed}
- Organellar genes quantified independently: {n_org_genes}
{"- MSTRG→organellar gene mapping: " + str(len(org_mapping)) + " genes mapped" if org_mapping is not None and not org_mapping.empty else ""}
"""
    with open(output_dir / "organellar_analysis_results.md", "w") as f:
        f.write(narrative)
    logger.info("Wrote organellar_analysis_results.md")


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    REPORT_PKG02.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {REPORT_PKG02}")

    # --- Load pipeline outputs ---
    # These paths may need adjustment based on actual pipeline output locations
    org_count_path = None
    org_de_path = None
    contam_path = None
    org_gff_paths = []
    org_fasta_paths = []

    # Search for outputs in pipeline results or reports directory
    search_dirs = [PIPELINE_RESULTS, MLIB / "reports"]
    for sd in search_dirs:
        if not sd.exists():
            continue
        for p in sd.rglob("*.csv"):
            if "organellar_gene_count_table" in p.name or "organellar_gene_count_table" in str(p):
                org_count_path = org_count_path or p
            if "deseq2_results" in p.name or "deseq2_results" in str(p):
                org_de_path = org_de_path or p
        for p in sd.rglob("*.tsv"):
            if "contamination_report" in p.name or "contamination" in str(p):
                contam_path = contam_path or p
        for p in sd.rglob("*.gff3"):
            org_gff_paths.append(p)
        for p in sd.rglob("*.fna"):
            org_fasta_paths.append(p)

    # --- Contamination report ---
    if contam_path:
        blast_df, contaminated_contigs = load_contamination_report(contam_path)
        shutil.copy(contam_path, REPORT_PKG02 / "contamination_report.tsv")
    else:
        logger.warning("Contamination report not found — assuming no contamination")
        blast_df = pd.DataFrame(columns=BLAST_HEADER)
        contaminated_contigs = set()

    # --- Gene ID mapping from report_package01 ---
    mapping_path = REPORT_PKG01 / "gene_id_mapping.csv"
    if not mapping_path.exists():
        # Try enrichment results
        mapping_path = RESULTS_BASE.parent / "enrichment" / "gene_id_mapping.csv"
    if mapping_path.exists():
        mapping_df = load_gene_id_mapping(mapping_path)
        organellar_mstrg_genes = identify_organellar_mstrg_genes(mapping_df, contaminated_contigs)
    else:
        logger.warning(f"Gene ID mapping not found at {mapping_path}")
        mapping_df = pd.DataFrame()
        organellar_mstrg_genes = set()

    # --- Build organellar gene mapping ---
    org_mapping = build_organellar_gene_mapping(blast_df, mapping_df, contaminated_contigs)
    if not org_mapping.empty:
        org_mapping.to_csv(REPORT_PKG02 / "organellar_gene_mapping.csv", index=False)

    # --- Clean nuclear tables ---
    nuc_count_path = REPORT_PKG01 / "gene_count_table.csv"
    nuc_de_path = REPORT_PKG01 / "gene_de_results.csv"
    n_removed = len(organellar_mstrg_genes)

    if nuc_count_path.exists() and nuc_de_path.exists():
        clean_nuclear_tables(nuc_count_path, nuc_de_path, organellar_mstrg_genes, REPORT_PKG02)
    else:
        logger.warning("Nuclear count/DE tables not found in report_package01")
        # Copy empty placeholders
        n_removed = 0

    # --- Copy organellar pipeline outputs ---
    org_counts = None
    org_de = None

    if org_count_path and org_count_path.exists():
        org_counts = pd.read_csv(org_count_path, index_col="gene_id")
        shutil.copy(org_count_path, REPORT_PKG02 / "organellar_gene_count_table.csv")
        logger.info(f"Organellar count table: {len(org_counts)} genes")
    else:
        logger.warning("Organellar count table not found")

    if org_de_path and org_de_path.exists():
        org_de = pd.read_csv(org_de_path, index_col="gene_id")
        shutil.copy(org_de_path, REPORT_PKG02 / "organellar_de_results.csv")
        logger.info(f"Organellar DE results: {len(org_de)} genes")
    else:
        logger.warning("Organellar DE results not found")

    # --- Copy GFF3 and FASTA files ---
    for gff_path in org_gff_paths:
        # Name by organelle if possible
        name = gff_path.stem
        if "chloroplast" in str(gff_path).lower() or "nc_023133" in str(gff_path).lower():
            name = "chloroplast"
        elif "mitochond" in str(gff_path).lower() or "mt483997" in str(gff_path).lower():
            name = "mitochondrion"
        shutil.copy(gff_path, REPORT_PKG02 / f"{name}.gff3")
        logger.info(f"Copied GFF3: {name}.gff3")

    for fasta_path in org_fasta_paths:
        name = fasta_path.stem
        if "chloroplast" in str(fasta_path).lower() or "nc_023133" in str(fasta_path).lower():
            name = "chloroplast"
        elif "mitochond" in str(fasta_path).lower() or "mt483997" in str(fasta_path).lower():
            name = "mitochondrion"
        shutil.copy(fasta_path, REPORT_PKG02 / f"{name}.fna")
        logger.info(f"Copied FASTA: {name}.fna")

    # --- Copy DAG images ---
    dag_dir = MLIB / "reports"
    for ext in ("svg", "png"):
        dag_file = dag_dir / f"organellar_workflow_dag.{ext}"
        if dag_file.exists():
            shutil.copy(dag_file, REPORT_PKG02 / dag_file.name)
            logger.info(f"Copied DAG: {dag_file.name}")

    # --- Compile report documents ---
    n_org_genes = len(org_counts) if org_counts is not None else 0
    compile_methods(REPORT_PKG02, n_org_genes, len(contaminated_contigs), n_removed)
    compile_results_narrative(
        REPORT_PKG02, org_counts, org_de, blast_df,
        contaminated_contigs, n_removed, org_mapping,
    )

    # --- Summary ---
    print("\n" + "=" * 60)
    print("report_package02 contents:")
    print("=" * 60)
    for p in sorted(REPORT_PKG02.iterdir()):
        size = p.stat().st_size
        if size > 1024 * 1024:
            size_str = f"{size / 1024 / 1024:.1f} MB"
        elif size > 1024:
            size_str = f"{size / 1024:.1f} KB"
        else:
            size_str = f"{size} B"
        print(f"  {p.name:45s} {size_str}")
    print(f"\nOutput: {REPORT_PKG02}")


if __name__ == "__main__":
    main()
