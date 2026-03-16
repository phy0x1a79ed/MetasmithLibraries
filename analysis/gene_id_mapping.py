"""Map MSTRG gene IDs (StringTie) to BRAKER3 gene IDs via genomic coordinate overlap."""

import re
import csv
import logging
from pathlib import Path
from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree

logger = logging.getLogger(__name__)


def parse_gtf_attribute(attr_str, key):
    """Extract value for a key from GTF attribute string."""
    match = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return match.group(1) if match else None


def parse_mstrg_genes(gtf_path):
    """Parse merged.gtf to get MSTRG gene coordinate spans.

    Returns dict: gene_id -> (contig, strand, start, end)
    """
    gene_spans = {}  # gene_id -> [contig, strand, min_start, max_end]

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            contig = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]

            gene_id = parse_gtf_attribute(attrs, 'gene_id')
            if not gene_id:
                continue

            if gene_id in gene_spans:
                gene_spans[gene_id][2] = min(gene_spans[gene_id][2], start)
                gene_spans[gene_id][3] = max(gene_spans[gene_id][3], end)
            else:
                gene_spans[gene_id] = [contig, strand, start, end]

    return {gid: tuple(v) for gid, v in gene_spans.items()}


def parse_braker3_genes(gff3_path):
    """Parse braker3.gff3 for gene-level coordinates.

    Returns dict: gene_id -> (contig, strand, start, end)
    """
    genes = {}

    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue

            contig = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]

            # Extract ID from GFF3 attributes
            match = re.search(r'ID=([^;]+)', attrs)
            if match:
                gene_id = match.group(1)
                genes[gene_id] = (contig, strand, start, end)

    return genes


def build_coordinate_mapping(mstrg_genes, braker3_genes, overlap_threshold=0.5):
    """Map MSTRG genes to BRAKER3 genes using coordinate overlap.

    Uses intervaltree for efficient overlap finding.
    Resolves ambiguities with reciprocal overlap fraction.
    """
    # Build interval trees per (contig, strand)
    trees = defaultdict(lambda: IntervalTree())
    for gene_id, (contig, strand, start, end) in braker3_genes.items():
        # intervaltree uses half-open intervals [start, end)
        trees[(contig, strand)].addi(start, end + 1, gene_id)

    mapping = {}  # mstrg_id -> (braker3_id, overlap_fraction)
    stats = {'total_mstrg': len(mstrg_genes), 'mapped': 0, 'unmapped': 0, 'ambiguous': 0}

    for mstrg_id, (contig, strand, m_start, m_end) in mstrg_genes.items():
        key = (contig, strand)
        if key not in trees:
            stats['unmapped'] += 1
            continue

        overlaps = trees[key].overlap(m_start, m_end + 1)
        if not overlaps:
            stats['unmapped'] += 1
            continue

        # Calculate reciprocal overlap for each candidate
        best_id = None
        best_score = 0

        for iv in overlaps:
            b_start, b_end, b_id = iv.begin, iv.end - 1, iv.data

            overlap_start = max(m_start, b_start)
            overlap_end = min(m_end, b_end)
            overlap_len = overlap_end - overlap_start + 1

            m_len = m_end - m_start + 1
            b_len = b_end - b_start + 1

            # Reciprocal overlap: min of both fractions
            recip_overlap = min(overlap_len / m_len, overlap_len / b_len)

            if recip_overlap > best_score:
                best_score = recip_overlap
                best_id = b_id

        if best_score >= overlap_threshold:
            mapping[mstrg_id] = (best_id, best_score)
            stats['mapped'] += 1
        elif best_id is not None:
            # Below threshold but has some overlap — still map but flag
            mapping[mstrg_id] = (best_id, best_score)
            stats['mapped'] += 1
            stats['ambiguous'] += 1
        else:
            stats['unmapped'] += 1

    return mapping, stats


def load_eggnog_annotations(eggnog_path):
    """Load eggNOG annotations indexed by query ID (g*.t*)."""
    rows = []
    with open(eggnog_path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#query'):
                # Header line
                headers = line.lstrip('#').strip().split('\t')
                continue
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= len(headers):
                rows.append(dict(zip(headers, fields)))

    df = pd.DataFrame(rows)
    df = df.set_index('query')
    return df


def parse_fasta_ids(fasta_path):
    """Parse FASTA headers and return set of sequence IDs."""
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].strip().split()[0])
    return ids


def build_full_mapping(report_dir, output_dir):
    """Build complete gene ID mapping and save results.

    Returns: (mapping_df, eggnog_df, stats)
    """
    report_dir = Path(report_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Parsing MSTRG genes from merged.gtf...")
    mstrg_genes = parse_mstrg_genes(report_dir / 'merged.gtf')
    logger.info(f"  Found {len(mstrg_genes)} MSTRG genes")

    logger.info("Parsing BRAKER3 genes from braker3.gff3...")
    braker3_genes = parse_braker3_genes(report_dir / 'braker3.gff3')
    logger.info(f"  Found {len(braker3_genes)} BRAKER3 genes")

    logger.info("Building coordinate mapping...")
    mapping, stats = build_coordinate_mapping(mstrg_genes, braker3_genes)
    logger.info(f"  Mapped: {stats['mapped']}/{stats['total_mstrg']} "
                f"({100*stats['mapped']/stats['total_mstrg']:.1f}%)")
    logger.info(f"  Unmapped: {stats['unmapped']}, Ambiguous: {stats['ambiguous']}")

    logger.info("Loading eggNOG annotations...")
    eggnog_df = load_eggnog_annotations(report_dir / 'eggnog_results.tsv')
    logger.info(f"  Loaded {len(eggnog_df)} annotations")

    # Load protein FASTA IDs (gffread output fed to eggNOG-mapper)
    protein_fasta_ids = set()
    protein_fasta_path = report_dir / 'proteins.faa'
    if protein_fasta_path.exists():
        protein_fasta_ids = parse_fasta_ids(protein_fasta_path)
        logger.info(f"  Loaded {len(protein_fasta_ids)} protein FASTA IDs")

    # Build mapping dataframe with coordinates
    rows = []
    for mstrg_id, (braker3_id, overlap_frac) in mapping.items():
        transcript_id = braker3_id + '.t1'
        has_annotation = transcript_id in eggnog_df.index
        in_fasta = transcript_id in protein_fasta_ids
        m_contig, m_strand, m_start, m_end = mstrg_genes[mstrg_id]
        b_contig, b_strand, b_start, b_end = braker3_genes[braker3_id]
        rows.append({
            'mstrg_id': mstrg_id,
            'braker3_gene_id': braker3_id,
            'braker3_transcript_id': transcript_id,
            'protein_fasta_id': transcript_id if in_fasta else '',
            'in_protein_fasta': in_fasta,
            'contig': m_contig,
            'strand': m_strand,
            'mstrg_start': m_start,
            'mstrg_end': m_end,
            'braker3_start': b_start,
            'braker3_end': b_end,
            'overlap_fraction': round(overlap_frac, 4),
            'has_eggnog_annotation': has_annotation,
        })

    mapping_df = pd.DataFrame(rows)

    # Add annotation columns for mapped genes
    annotated = mapping_df[mapping_df['has_eggnog_annotation']].copy()
    if len(annotated) > 0:
        ann_cols = ['GOs', 'KEGG_ko', 'KEGG_Pathway', 'COG_category', 'Description', 'PFAMs']
        ann_data = eggnog_df.loc[annotated['braker3_transcript_id'],
                                  [c for c in ann_cols if c in eggnog_df.columns]]
        ann_data.index = annotated.index
        for col in ann_data.columns:
            mapping_df.loc[annotated.index, col] = ann_data[col]

    # Save
    out_path = output_dir / 'gene_id_mapping.csv'
    mapping_df.to_csv(out_path, index=False)
    logger.info(f"Saved mapping to {out_path}")

    # Stats summary
    n_with_ann = mapping_df['has_eggnog_annotation'].sum()
    stats['with_annotation'] = int(n_with_ann)
    stats['braker3_total'] = len(braker3_genes)
    stats['eggnog_total'] = len(eggnog_df)

    return mapping_df, eggnog_df, stats


def export_clean_mapping(mapping_df, output_path):
    """Export a clean mapping table with coordinates and description only.

    Drops heavy annotation columns (GOs, KEGG, COG, PFAMs) — keeps just
    the mapping, coordinates, and a human-readable description.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    clean_cols = [
        'mstrg_id', 'braker3_gene_id', 'braker3_transcript_id',
        'protein_fasta_id', 'in_protein_fasta',
        'contig', 'strand', 'mstrg_start', 'mstrg_end',
        'braker3_start', 'braker3_end', 'overlap_fraction',
        'has_eggnog_annotation', 'eggnog_description',
    ]

    clean_df = mapping_df.copy()
    # Rename Description -> eggnog_description if present
    if 'Description' in clean_df.columns:
        clean_df = clean_df.rename(columns={'Description': 'eggnog_description'})

    # Keep only columns that exist
    clean_df = clean_df[[c for c in clean_cols if c in clean_df.columns]]

    clean_df.to_csv(output_path, index=False)
    logger.info(f"Exported clean mapping ({len(clean_df)} rows) to {output_path}")
    return clean_df


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    build_full_mapping(
        'results/report_package01',
        'results/enrichment'
    )
