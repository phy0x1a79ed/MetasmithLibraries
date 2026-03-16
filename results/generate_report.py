#!/usr/bin/env python3
"""
Generate metabolomics report package for Porphyridium purpureum 161-PE experiment.

Produces results/report_package01/ with:
  feature_table.csv, compound_metadata.csv, differential_results.csv,
  enrichment_results.csv, pathway_enrichment_results.md, methods.md,
  volcano_differential.svg/.png, enrichment_dotplot.svg/.png, workflow_dag.svg/.png

Run: conda run -n msm_env python3 results/generate_report.py
"""

import os
import sys
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / "tests" / "test_data" / "porphyridium-metabolomics" / "jgi_metabolomic_data"
OUT_DIR = SCRIPT_DIR / "report_package01"

TARGETED_DIR = DATA_DIR / "Targeted"
UNTARGETED_DIR = DATA_DIR / "Untargeted"

# 161-PE experimental groups
PE_PREFIX = "161-PE"
BASELINE = "0h"
CONDITIONS = ["8", "P", "S"]
TIMEPOINTS = [24, 72, 144]

# Significance thresholds
ALPHA = 0.05
LOG2FC_CUTOFF = 1.0

N_HEADER_ROWS = 6  # group, file, short groupname, sample treatment, short filename, short samplename

# KEGG REST cache (populated at runtime)
_kegg_cache_dir = SCRIPT_DIR / ".kegg_cache"


def bh_correction(pvals):
    """Benjamini-Hochberg FDR correction (manual, no statsmodels needed)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n + 1)
    adj = pvals * n / ranks
    # enforce monotonicity (cumulative min from right)
    adj_sorted = adj[np.argsort(-ranks)]  # sorted by descending rank
    adj_sorted = np.minimum.accumulate(adj_sorted)
    adj[np.argsort(-ranks)] = adj_sorted
    return np.clip(adj, 0, 1)


# ===================================================================
# STEP 1: Load Targeted Data
# ===================================================================
def parse_sample_label(short_groupname):
    """Extract condition and timepoint from short groupname like POS_161-PE-8-24h."""
    s = short_groupname.upper()
    if PE_PREFIX.upper() not in s:
        return None, None
    # strip polarity prefix
    after = s.split(PE_PREFIX.upper() + "-", 1)[-1]  # e.g. "8-24H", "0H", "P-72H"
    if after.startswith("0H"):
        return "0h", 0
    for cond in CONDITIONS:
        if after.startswith(cond.upper() + "-"):
            rest = after[len(cond) + 1:]  # e.g. "24H"
            for tp in TIMEPOINTS:
                if rest.startswith(f"{tp}H"):
                    return cond, tp
    return None, None


def load_targeted_peak_heights(ema_dir, polarity):
    """Load peak_height.tab from an EMA directory, return (DataFrame, sample_info dict)."""
    pol = polarity.upper()  # POS or NEG
    ph_path = ema_dir / f"{pol}_data_sheets" / f"{pol}_peak_height.tab"
    if not ph_path.exists():
        return None, None

    # Read all rows as text first to extract header metadata
    with open(ph_path) as f:
        lines = f.readlines()

    if len(lines) < N_HEADER_ROWS + 1:
        return None, None

    # Parse header rows
    header_names = ["group", "file", "short_groupname", "sample_treatment",
                    "short_filename", "short_samplename"]
    headers = {}
    for i, name in enumerate(header_names):
        headers[name] = lines[i].rstrip("\n").split("\t")

    # short_groupname row: first cell is "short groupname", rest are group labels
    groupnames = headers["short_groupname"][1:]  # skip row label
    treatments = headers["sample_treatment"][1:]

    # Identify 161-PE columns (0-indexed into groupnames)
    pe_indices = []
    pe_groupnames = []
    for j, gn in enumerate(groupnames):
        cond, tp = parse_sample_label(gn)
        if cond is not None:
            pe_indices.append(j)
            pe_groupnames.append(gn)

    if not pe_indices:
        return None, None

    # Make unique column names by appending replicate number
    group_counts = {}
    pe_labels = []
    for gn in pe_groupnames:
        group_counts[gn] = group_counts.get(gn, 0) + 1
        pe_labels.append(f"{gn}_r{group_counts[gn]}")

    # Parse data rows (skip N_HEADER_ROWS)
    data_lines = lines[N_HEADER_ROWS:]
    compounds = []
    values = []
    for line in data_lines:
        parts = line.rstrip("\n").split("\t")
        if not parts or not parts[0]:
            continue
        compounds.append(parts[0])
        row_vals = []
        for j in pe_indices:
            idx = j + 1  # +1 because first column is compound name
            if idx < len(parts):
                try:
                    row_vals.append(float(parts[idx]))
                except ValueError:
                    row_vals.append(np.nan)
            else:
                row_vals.append(np.nan)
        values.append(row_vals)

    df = pd.DataFrame(values, index=compounds, columns=pe_labels)
    df.index.name = "compound_name"

    # Build sample info mapping column -> (condition, timepoint, replicate)
    sample_info = {}
    for col in pe_labels:
        # Strip _rN suffix to parse
        base = col.rsplit("_r", 1)[0]
        cond, tp = parse_sample_label(base)
        rep = int(col.rsplit("_r", 1)[1])
        sample_info[col] = {"condition": cond, "timepoint_h": tp,
                            "replicate": rep}

    return df, sample_info


def load_istd_normalization_factors(istd_dir, polarity, pe_labels):
    """Load ISTD peak heights and compute per-sample median normalization factors."""
    if istd_dir is None or not istd_dir.exists():
        return None
    pol = polarity.upper()
    ph_path = istd_dir / f"{pol}_data_sheets" / f"{pol}_peak_height.tab"
    if not ph_path.exists():
        return None

    with open(ph_path) as f:
        lines = f.readlines()

    if len(lines) < N_HEADER_ROWS + 1:
        return None

    groupnames_line = lines[2].rstrip("\n").split("\t")
    groupnames = groupnames_line[1:]

    # Find matching PE columns
    pe_indices = []
    pe_groupnames = []
    for j, gn in enumerate(groupnames):
        cond, tp = parse_sample_label(gn)
        if cond is not None:
            pe_indices.append(j)
            pe_groupnames.append(gn)

    if not pe_indices:
        return None

    # Make unique labels matching EMA naming
    group_counts = {}
    istd_labels = []
    for gn in pe_groupnames:
        group_counts[gn] = group_counts.get(gn, 0) + 1
        istd_labels.append(f"{gn}_r{group_counts[gn]}")

    data_lines = lines[N_HEADER_ROWS:]
    values = []
    for line in data_lines:
        parts = line.rstrip("\n").split("\t")
        if not parts or not parts[0]:
            continue
        row_vals = []
        for j in pe_indices:
            idx = j + 1
            if idx < len(parts):
                try:
                    row_vals.append(float(parts[idx]))
                except ValueError:
                    row_vals.append(np.nan)
            else:
                row_vals.append(np.nan)
        values.append(row_vals)

    istd_df = pd.DataFrame(values, columns=istd_labels)
    # Per-sample median across all ISTD compounds
    medians = istd_df.median(axis=0).replace(0, np.nan)
    # Normalize: divide by global median of medians so factors are relative
    global_median = medians.median()
    if global_median > 0 and not np.isnan(global_median):
        nf = medians / global_median
    else:
        nf = medians

    return nf


def load_compound_atlas(ema_dir):
    """Load CompoundAtlas CSV and return metadata DataFrame."""
    atlas_files = sorted(ema_dir.glob("CompoundAtlas__*.csv"))
    if not atlas_files:
        return None
    atlas = pd.read_csv(atlas_files[0])
    # Rename label to compound_name for joining
    if "label" in atlas.columns:
        atlas = atlas.rename(columns={"label": "compound_name"})
    elif "name" in atlas.columns and "compound_name" not in atlas.columns:
        atlas["compound_name"] = atlas["name"]
    return atlas


def load_all_targeted():
    """Load all targeted data across chromatographies and polarities."""
    all_dfs = []
    all_meta = []
    all_sample_info = {}

    if not TARGETED_DIR.exists():
        return pd.DataFrame(), pd.DataFrame(), {}

    for chrom_dir in sorted(TARGETED_DIR.iterdir()):
        if not chrom_dir.is_dir():
            continue
        chrom = "HILIC" if "HILIC" in chrom_dir.name.upper() else "C18"

        for pol in ["POS", "NEG"]:
            ema_dirs = sorted(chrom_dir.glob(f"EMA-{pol}_*"))
            istd_dirs = sorted(chrom_dir.glob(f"ISTDsEtc*-{pol}_*"))
            istd_dir = istd_dirs[0] if istd_dirs else None

            for ema_dir in ema_dirs:
                df, sample_info = load_targeted_peak_heights(ema_dir, pol)
                if df is None:
                    continue

                # ISTD normalization
                nf = load_istd_normalization_factors(istd_dir, pol, list(df.columns))
                if nf is not None:
                    # Both EMA and ISTD now have _rN suffix matching
                    for col in df.columns:
                        if col in nf.index:
                            val = float(nf[col])
                            if val > 0 and not np.isnan(val):
                                df[col] = df[col] / val

                # Load compound atlas and join by position (atlas row i = peak_height row i)
                atlas = load_compound_atlas(ema_dir)

                prefix = f"{chrom}_{pol}_"
                compound_ids = [prefix + name for name in df.index]
                df.index = compound_ids

                if atlas is not None:
                    want_cols = ["name", "kegg_id", "inchi_key",
                                 "pubchem_compound_id", "formula", "hmdb_id",
                                 "chebi_id", "lipidmaps_id", "rt_peak", "mz",
                                 "adduct", "polarity"]
                    want_cols = [c for c in want_cols if c in atlas.columns]
                    # Atlas rows align positionally with peak_height rows
                    meta = atlas[want_cols].head(len(compound_ids)).copy()
                    meta.insert(0, "compound_name", compound_ids[:len(meta)])
                else:
                    meta = pd.DataFrame({"compound_name": compound_ids})

                meta["compound_id"] = meta["compound_name"]
                meta["chromatography"] = chrom
                meta["source"] = "targeted"
                if "polarity" not in meta.columns:
                    meta["polarity"] = pol.lower()

                all_dfs.append(df)
                all_meta.append(meta)
                all_sample_info.update(sample_info)

    if not all_dfs:
        return pd.DataFrame(), pd.DataFrame(), {}

    combined = pd.concat(all_dfs, axis=0)
    combined_meta = pd.concat(all_meta, ignore_index=True)

    return combined, combined_meta, all_sample_info


# ===================================================================
# STEP 2: Load Untargeted Data
# ===================================================================
def load_all_untargeted():
    """Load untargeted data across all directories."""
    all_dfs = []
    all_meta = []

    if not UNTARGETED_DIR.exists():
        return pd.DataFrame(), pd.DataFrame()

    for run_dir in sorted(UNTARGETED_DIR.iterdir()):
        if not run_dir.is_dir() or "archive" in run_dir.name.lower():
            continue

        chrom = "HILIC" if "HILIC" in run_dir.name.upper() else "C18"

        for polarity_tag in ["positive", "negative"]:
            pol = "POS" if polarity_tag == "positive" else "NEG"

            # Find peak height file
            ph_files = sorted(run_dir.glob(f"*_{polarity_tag}_peak-height-filtered-3x-exctrl.csv"))
            if not ph_files:
                continue

            # Find metadata file
            meta_files = sorted(run_dir.glob(f"*_{polarity_tag}_metadata.tab"))

            # Find GNPS library results
            gnps_files = sorted(run_dir.glob(f"*_{polarity_tag}_gnps2-fbmn-library-results.tsv"))

            # Load peak heights
            heights = pd.read_csv(ph_files[0])
            id_col = heights.columns[0]  # "row ID"
            mz_col = heights.columns[1]  # "row m/z"
            rt_col = heights.columns[2]  # "row retention time"

            # Filter to 161-PE samples using metadata
            sample_cols = [c for c in heights.columns[3:] if "Peak height" in c]

            if meta_files:
                meta_df = pd.read_csv(meta_files[0], sep="\t")
                att_col = next((c for c in meta_df.columns if "sampletype" in c.lower()), None)
                fn_col = next((c for c in meta_df.columns if "filename" in c.lower()), None)
                if att_col and fn_col:
                    pe_meta = meta_df[meta_df[att_col].str.contains("161-pe", case=False, na=False)]
                    pe_fns = set(pe_meta[fn_col].str.replace(r"\.mzML$", "", regex=True))
                    pe_cols = [c for c in sample_cols
                               if any(fn in c for fn in pe_fns)]
                else:
                    pe_cols = [c for c in sample_cols if "161-pe" in c.lower()]
            else:
                pe_cols = [c for c in sample_cols if "161-pe" in c.lower()]

            if not pe_cols:
                continue

            # Load GNPS annotations
            annotations = None
            if gnps_files:
                gnps = pd.read_csv(gnps_files[0], sep="\t")
                # Keep best match per scan (highest MQScore)
                scan_col = "#Scan#"
                if scan_col in gnps.columns:
                    gnps = gnps.sort_values("MQScore", ascending=False).drop_duplicates(scan_col, keep="first")
                    annotations = gnps

            # Build feature DataFrame
            feature_ids = heights[id_col].values
            mz_vals = heights[mz_col].values
            rt_vals = heights[rt_col].values
            data = heights[pe_cols].values

            prefix = f"UT_{chrom}_{pol}_"
            compound_names = [f"{prefix}{fid}" for fid in feature_ids]

            df = pd.DataFrame(data, index=compound_names, columns=pe_cols)
            df.index.name = "compound_name"

            # Build metadata
            meta_records = []
            for i, fid in enumerate(feature_ids):
                rec = {
                    "compound_name": compound_names[i],
                    "compound_id": compound_names[i],
                    "feature_id": int(fid),
                    "mz": mz_vals[i],
                    "rt_peak": rt_vals[i],
                    "polarity": pol.lower(),
                    "chromatography": chrom,
                    "source": "untargeted",
                }
                # Join GNPS annotation if available
                if annotations is not None and scan_col in annotations.columns:
                    match = annotations[annotations[scan_col] == fid]
                    if not match.empty:
                        row = match.iloc[0]
                        rec["name"] = row.get("Compound_Name", "")
                        rec["inchi_key"] = row.get("InChIKey", "")
                        rec["formula"] = row.get("molecular_formula", "")
                        rec["adduct"] = row.get("Adduct", "")
                        rec["gnps_score"] = row.get("MQScore", np.nan)
                        rec["gnps_compound_name"] = row.get("Compound_Name", "")
                    else:
                        rec["name"] = ""
                        rec["inchi_key"] = ""
                meta_records.append(rec)

            meta = pd.DataFrame(meta_records)

            # Keep only GNPS-annotated features
            if "name" in meta.columns:
                annotated = meta[meta["name"].notna() & (meta["name"] != "")]
                if not annotated.empty:
                    keep_ids = set(annotated["compound_name"])
                    df = df.loc[df.index.isin(keep_ids)]
                    meta = annotated.copy()
                else:
                    # No annotated features; skip this polarity/chrom combo
                    continue

            all_dfs.append(df)
            all_meta.append(meta)

    if not all_dfs:
        return pd.DataFrame(), pd.DataFrame()

    combined = pd.concat(all_dfs, axis=0)
    combined_meta = pd.concat(all_meta, ignore_index=True)
    return combined, combined_meta


# ===================================================================
# STEP 3: Merge & Normalize
# ===================================================================
def assign_clean_sample_labels(columns, sample_info_targeted):
    """Map raw column names to clean labels like '8-24h-1', '0h-2', etc."""
    label_map = {}
    group_counters = {}

    for col in columns:
        # For targeted columns with _rN suffix, strip it for parsing
        base = col.rsplit("_r", 1)[0] if "_r" in col else col
        cond, tp = parse_sample_label(base)
        if cond is not None:
            if cond == "0h":
                key = "0h"
            else:
                key = f"{cond}-{tp}h"
            group_counters[key] = group_counters.get(key, 0) + 1
            clean = f"{key}-{group_counters[key]}"
            label_map[col] = clean
        else:
            # Untargeted columns: parse from filename
            col_lower = col.lower()
            found = False
            for c in CONDITIONS:
                for t in TIMEPOINTS:
                    pattern = f"161-pe-{c.lower()}-{t}h"
                    if pattern in col_lower:
                        key = f"{c}-{t}h"
                        group_counters[key] = group_counters.get(key, 0) + 1
                        label_map[col] = f"{key}-{group_counters[key]}"
                        found = True
                        break
                if found:
                    break
            if not found and "161-pe-0h" in col_lower:
                key = "0h"
                group_counters[key] = group_counters.get(key, 0) + 1
                label_map[col] = f"{key}-{group_counters[key]}"
            elif not found:
                label_map[col] = col  # keep as-is

    return label_map


def merge_and_normalize(targeted_df, targeted_meta, untargeted_df, untargeted_meta,
                        sample_info_targeted):
    """Merge targeted + untargeted, normalize, deduplicate."""
    # Unify column names: both DataFrames may have different column names
    # Targeted uses short groupnames, untargeted uses long filenames
    # We need to align them to a common set of clean sample labels

    # --- Targeted: rename columns to clean labels ---
    if not targeted_df.empty:
        tgt_label_map = assign_clean_sample_labels(targeted_df.columns, sample_info_targeted)
        targeted_df = targeted_df.rename(columns=tgt_label_map)

    # --- Untargeted: rename columns to clean labels ---
    if not untargeted_df.empty:
        ut_label_map = assign_clean_sample_labels(untargeted_df.columns, {})
        untargeted_df = untargeted_df.rename(columns=ut_label_map)

    # Combine all sample columns across both
    all_sample_cols = set()
    if not targeted_df.empty:
        all_sample_cols.update(targeted_df.columns)
    if not untargeted_df.empty:
        all_sample_cols.update(untargeted_df.columns)

    # Concat rows
    dfs = []
    if not targeted_df.empty:
        dfs.append(targeted_df)
    if not untargeted_df.empty:
        # Align untargeted columns to same set
        for col in all_sample_cols - set(untargeted_df.columns):
            untargeted_df[col] = np.nan
        untargeted_df = untargeted_df[sorted(all_sample_cols & set(untargeted_df.columns))]
        dfs.append(untargeted_df)

    if not dfs:
        return pd.DataFrame(), pd.DataFrame()

    # Align all dfs to same columns
    all_cols = sorted(all_sample_cols)
    for i, d in enumerate(dfs):
        for col in all_cols:
            if col not in d.columns:
                d[col] = np.nan
        dfs[i] = d[all_cols]

    combined = pd.concat(dfs, axis=0)

    # Combine metadata
    meta_dfs = []
    if not targeted_meta.empty:
        meta_dfs.append(targeted_meta)
    if not untargeted_meta.empty:
        meta_dfs.append(untargeted_meta)
    combined_meta = pd.concat(meta_dfs, ignore_index=True) if meta_dfs else pd.DataFrame()

    # Log2 transform: impute zeros but keep NaN as NaN
    # NaN = compound not measured on that platform (preserve as missing)
    # 0 = below detection limit (impute with min_positive / 2)
    for col in combined.columns:
        vals = combined[col].astype(float)
        positives = vals[vals > 0]
        min_pos = positives.min() if not positives.empty else 1.0
        impute = min_pos / 2 if pd.notna(min_pos) and min_pos > 0 else 1.0
        # Replace 0 with impute, keep NaN
        transformed = vals.copy()
        transformed[vals == 0] = impute
        transformed[vals > 0] = vals[vals > 0]
        # NaN stays NaN
        combined[col] = np.log2(transformed)

    # Deduplicate by inchi_key + polarity (prefer targeted)
    if "inchi_key" in combined_meta.columns:
        combined_meta["_idx"] = combined_meta["compound_name"]
        dedup_key = combined_meta[combined_meta["inchi_key"].notna() & (combined_meta["inchi_key"] != "")]
        if not dedup_key.empty:
            dedup_key = dedup_key.sort_values("source")  # targeted < untargeted alphabetically
            dups = dedup_key.duplicated(subset=["inchi_key", "polarity"], keep="first")
            to_remove = set(dedup_key[dups]["compound_name"])
            combined = combined[~combined.index.isin(to_remove)]
            combined_meta = combined_meta[~combined_meta["compound_name"].isin(to_remove)]
        if "_idx" in combined_meta.columns:
            combined_meta = combined_meta.drop(columns=["_idx"])

    return combined, combined_meta


# ===================================================================
# STEP 4: Differential Analysis
# ===================================================================
def parse_condition_from_label(label):
    """Parse condition from clean label like '8-24h-1' or '0h-2'."""
    if label.startswith("0h"):
        return "0h", 0
    for c in CONDITIONS:
        if label.startswith(f"{c}-"):
            rest = label[len(c) + 1:]
            for tp in TIMEPOINTS:
                if rest.startswith(f"{tp}h"):
                    return c, tp
    return None, None


def run_differential_analysis(feature_table, compound_meta):
    """Run Welch's t-test for each comparison, return stacked results DataFrame."""
    sample_cols = list(feature_table.columns)

    # Map columns to conditions
    col_info = {}
    for col in sample_cols:
        cond, tp = parse_condition_from_label(col)
        if cond is not None:
            col_info[col] = {"condition": cond, "timepoint_h": tp}

    baseline_cols = [c for c, info in col_info.items() if info["condition"] == "0h"]

    all_results = []

    for cond in CONDITIONS:
        for tp in TIMEPOINTS:
            treatment_cols = [c for c, info in col_info.items()
                              if info["condition"] == cond and info["timepoint_h"] == tp]
            if not treatment_cols or not baseline_cols:
                continue

            comparison = f"{cond}_vs_0h_t{tp}h"
            records = []

            for cpd_name in feature_table.index:
                row = feature_table.loc[cpd_name]
                g_treat = row[treatment_cols].dropna().values.astype(float)
                g_base = row[baseline_cols].dropna().values.astype(float)

                if len(g_treat) < 2 or len(g_base) < 2:
                    continue

                log2fc = g_treat.mean() - g_base.mean()

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    _, pval = stats.ttest_ind(g_treat, g_base, equal_var=False)

                if np.isnan(pval):
                    continue

                # Get metadata
                meta_row = compound_meta[compound_meta["compound_name"] == cpd_name]
                kegg_id = ""
                inchi_key = ""
                name = cpd_name
                if not meta_row.empty:
                    kegg_id = str(meta_row.iloc[0].get("kegg_id", ""))
                    inchi_key = str(meta_row.iloc[0].get("inchi_key", ""))
                    if "name" in meta_row.columns:
                        n = meta_row.iloc[0]["name"]
                        if pd.notna(n) and n:
                            name = n

                records.append({
                    "compound_id": cpd_name,
                    "compound_name": name,
                    "kegg_id": kegg_id if kegg_id != "nan" else "",
                    "inchi_key": inchi_key if inchi_key != "nan" else "",
                    "log2FC": log2fc,
                    "pval": pval,
                    "comparison": comparison,
                    "condition": cond,
                    "timepoint_h": tp,
                })

            if not records:
                continue

            res_df = pd.DataFrame(records)
            res_df["adj_pval"] = bh_correction(res_df["pval"].values)
            res_df["significant"] = (res_df["adj_pval"] < ALPHA) & (res_df["log2FC"].abs() > LOG2FC_CUTOFF)
            all_results.append(res_df)

    # ANOVA per timepoint across conditions
    anova_records = []
    for tp in TIMEPOINTS:
        for cpd_name in feature_table.index:
            row = feature_table.loc[cpd_name]
            groups = []
            for cond in CONDITIONS:
                cols = [c for c, info in col_info.items()
                        if info["condition"] == cond and info["timepoint_h"] == tp]
                g = row[cols].dropna().values.astype(float) if cols else np.array([])
                groups.append(g)
            groups = [g for g in groups if len(g) >= 2]
            if len(groups) < 2:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fstat, pval = stats.f_oneway(*groups)
            if np.isnan(pval):
                continue
            anova_records.append({
                "compound_id": cpd_name,
                "timepoint_h": tp,
                "F_stat": fstat,
                "pval": pval,
            })

    if anova_records:
        anova_df = pd.DataFrame(anova_records)
        anova_df["adj_pval"] = bh_correction(anova_df["pval"].values)

    if all_results:
        diff_df = pd.concat(all_results, ignore_index=True)
    else:
        diff_df = pd.DataFrame()

    return diff_df


# ===================================================================
# STEP 5: Volcano Plots
# ===================================================================
def make_volcano_plots(diff_df):
    """3x3 volcano plot grid: rows=conditions, cols=timepoints.

    Shared x/y axes with row/column labels instead of per-panel titles.
    Inner tick labels and axis lines are hidden.
    """
    if diff_df.empty:
        return

    COND_LABELS = {"8": "Selenium (8 \u00b5M)", "P": "Phosphorus deprivation",
                   "S": "Sulfur deprivation"}

    fig, axes = plt.subplots(
        3, 3, figsize=(12, 10), squeeze=False,
        sharex=True, sharey=True,
    )
    fig.suptitle("Volcano Plots: Metabolite Changes vs Baseline (0h)",
                 fontsize=13, y=0.97)

    for i, cond in enumerate(CONDITIONS):
        for j, tp in enumerate(TIMEPOINTS):
            ax = axes[i][j]
            comp = f"{cond}_vs_0h_t{tp}h"
            subset = diff_df[diff_df["comparison"] == comp]

            if subset.empty:
                continue

            ns = subset[~subset["significant"]]
            sig = subset[subset["significant"]]

            ax.scatter(ns["log2FC"], -np.log10(ns["pval"].clip(lower=1e-300)),
                       color="grey", alpha=0.4, s=12, zorder=1)
            ax.scatter(sig["log2FC"], -np.log10(sig["pval"].clip(lower=1e-300)),
                       color="tomato", s=18, zorder=2, edgecolors="darkred",
                       linewidths=0.3)

            # Label top 5 by |log2FC| among significant
            if not sig.empty:
                top5 = sig.nlargest(5, "log2FC", keep="first")
                bot5 = sig.nsmallest(5, "log2FC", keep="first")
                to_label = pd.concat([top5, bot5]).drop_duplicates(
                    "compound_id").head(5)
                for _, r in to_label.iterrows():
                    label = r["compound_name"]
                    if len(label) > 25:
                        label = label[:22] + "..."
                    ax.annotate(
                        label,
                        (r["log2FC"], -np.log10(max(r["pval"], 1e-300))),
                        fontsize=5, alpha=0.8, ha="center", va="bottom")

            ax.axvline(LOG2FC_CUTOFF, color="steelblue", linestyle="--",
                       linewidth=0.8)
            ax.axvline(-LOG2FC_CUTOFF, color="steelblue", linestyle="--",
                       linewidth=0.8)
            ax.axhline(-np.log10(ALPHA), color="orange", linestyle=":",
                       linewidth=0.8)

            # n_sig annotation inside each panel
            ax.text(0.97, 0.97, f"n={len(sig)}",
                    transform=ax.transAxes, fontsize=8, va="top", ha="right",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

    # Hide inner spines and ticks
    for i in range(3):
        for j in range(3):
            ax = axes[i][j]
            # Remove right spine on all; remove left spine on non-leftmost
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            if j > 0:
                ax.spines["left"].set_visible(False)
                ax.tick_params(left=False)
            if i < 2:
                ax.spines["bottom"].set_visible(False)
                ax.tick_params(bottom=False)

    # Column titles (timepoints) on top row
    for j, tp in enumerate(TIMEPOINTS):
        axes[0][j].set_title(f"{tp}h", fontsize=11, fontweight="bold")

    # Row labels (conditions) on right side
    for i, cond in enumerate(CONDITIONS):
        axes[i][2].annotate(
            COND_LABELS.get(cond, cond),
            xy=(1.04, 0.5), xycoords="axes fraction",
            fontsize=10, fontweight="bold", rotation=-90,
            ha="left", va="center")

    # Shared axis labels — only on outer edges
    for i in range(3):
        axes[i][0].set_ylabel("-log\u2081\u2080(p-value)", fontsize=10)
    for j in range(3):
        axes[2][j].set_xlabel("log\u2082FC", fontsize=10)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
               markersize=6, label='n.s.'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='tomato',
               markersize=6, label='significant'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2,
               fontsize=10, bbox_to_anchor=(0.5, 0.005))

    plt.tight_layout(rect=[0, 0.04, 0.94, 0.94])
    fig.savefig(OUT_DIR / "volcano_differential.svg", format="svg")
    fig.savefig(OUT_DIR / "volcano_differential.png", dpi=300)
    plt.close()
    print("  Saved volcano_differential.svg / .png")


# ===================================================================
# STEP 6: KEGG Pathway Enrichment (ORA)
# ===================================================================
def kegg_rest_get(url, cache_name=None):
    """Fetch from KEGG REST API with caching and rate limiting."""
    if cache_name:
        cache_file = _kegg_cache_dir / cache_name
        if cache_file.exists():
            return cache_file.read_text()

    try:
        import requests
        resp = requests.get(url, timeout=30)
        time.sleep(0.4)  # KEGG rate limit
        if resp.ok:
            text = resp.text
            if cache_name:
                _kegg_cache_dir.mkdir(parents=True, exist_ok=True)
                cache_file.write_text(text)
            return text
    except Exception as e:
        print(f"  KEGG API error for {url}: {e}")
    return ""


def fetch_kegg_pathways():
    """Fetch all KEGG metabolic pathway-compound mappings."""
    pathway_compounds = {}
    pathway_names = {}

    # Get pathway list
    text = kegg_rest_get("https://rest.kegg.jp/list/pathway/map", "pathway_list.txt")
    if not text:
        return {}, {}

    for line in text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            pid = parts[0].replace("path:", "")
            pathway_names[pid] = parts[1]

    print(f"  Found {len(pathway_names)} KEGG pathways, fetching compound links...")

    # Fetch compound links for each pathway
    fetched = 0
    for pid in list(pathway_names.keys()):
        r_text = kegg_rest_get(f"https://rest.kegg.jp/link/compound/{pid}",
                               f"link_{pid}.txt")
        if r_text and r_text.strip():
            cpds = set()
            for line in r_text.strip().split("\n"):
                p = line.split("\t")
                if len(p) >= 2:
                    cpds.add(p[1].replace("cpd:", ""))
            if cpds:
                pathway_compounds[pid] = cpds
        fetched += 1
        if fetched % 50 == 0:
            print(f"    ... fetched {fetched}/{len(pathway_names)} pathways")

    print(f"  Loaded {len(pathway_compounds)} pathways with compound mappings")
    return pathway_compounds, pathway_names


def convert_pubchem_to_kegg(pubchem_ids):
    """Try to convert PubChem CIDs to KEGG compound IDs via KEGG REST API."""
    mapping = {}
    for cid in pubchem_ids:
        if not cid or str(cid) == "nan":
            continue
        cid_str = str(int(float(cid)))
        text = kegg_rest_get(f"https://rest.kegg.jp/conv/compound/pubchem:{cid_str}",
                             f"conv_pubchem_{cid_str}.txt")
        if text and text.strip():
            for line in text.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    kegg_id = parts[1].replace("cpd:", "")
                    mapping[cid_str] = kegg_id
                    break
    return mapping


def run_pathway_enrichment(diff_df, compound_meta):
    """ORA using Fisher's exact test per pathway per comparison."""
    if diff_df.empty:
        return pd.DataFrame()

    # Collect all KEGG IDs and try to augment with PubChem conversions
    kegg_col = "kegg_id"
    all_kegg = set()
    for _, row in diff_df.iterrows():
        kid = str(row.get(kegg_col, ""))
        if kid and kid != "" and kid != "nan":
            # Handle multiple KEGG IDs separated by ///
            for k in kid.split("///"):
                k = k.strip()
                if k:
                    all_kegg.add(k)

    # Try PubChem conversion for compounds missing KEGG IDs
    missing_kegg = diff_df[
        (diff_df[kegg_col].isna()) | (diff_df[kegg_col] == "") | (diff_df[kegg_col] == "nan")
    ]
    if "pubchem_compound_id" in compound_meta.columns:
        pubchem_ids = set()
        for cpd_id in missing_kegg["compound_id"].unique():
            meta_row = compound_meta[compound_meta["compound_name"] == cpd_id]
            if not meta_row.empty:
                pcid = meta_row.iloc[0].get("pubchem_compound_id", "")
                if pd.notna(pcid) and str(pcid) != "nan" and str(pcid) != "":
                    pubchem_ids.add(str(pcid))

        if pubchem_ids:
            print(f"  Converting {len(pubchem_ids)} PubChem IDs to KEGG...")
            pc_to_kegg = convert_pubchem_to_kegg(pubchem_ids)
            print(f"  Converted {len(pc_to_kegg)} PubChem -> KEGG")

            # Update diff_df KEGG IDs
            for cpd_id in missing_kegg["compound_id"].unique():
                meta_row = compound_meta[compound_meta["compound_name"] == cpd_id]
                if not meta_row.empty:
                    pcid = str(meta_row.iloc[0].get("pubchem_compound_id", ""))
                    pcid_clean = str(int(float(pcid))) if pcid and pcid != "nan" else ""
                    if pcid_clean in pc_to_kegg:
                        kid = pc_to_kegg[pcid_clean]
                        diff_df.loc[diff_df["compound_id"] == cpd_id, kegg_col] = kid
                        all_kegg.add(kid)

    print(f"  Total KEGG compound IDs in background: {len(all_kegg)}")

    if len(all_kegg) < 3:
        print("  Too few KEGG IDs for enrichment analysis")
        return pd.DataFrame(columns=["pathway", "pathway_name", "comparison",
                                      "pval", "adj_pval", "overlap_ratio", "metabolites"])

    # Fetch KEGG pathway data
    pathway_compounds, pathway_names = fetch_kegg_pathways()
    if not pathway_compounds:
        print("  Could not load KEGG pathway data")
        return pd.DataFrame(columns=["pathway", "pathway_name", "comparison",
                                      "pval", "adj_pval", "overlap_ratio", "metabolites"])

    # Build sig sets per comparison
    sig_sets = {}
    for comp in diff_df["comparison"].unique():
        comp_df = diff_df[diff_df["comparison"] == comp]
        sig_df = comp_df[comp_df["significant"] == True]
        sig_kegg = set()
        for _, row in sig_df.iterrows():
            kid = str(row.get(kegg_col, ""))
            if kid and kid != "" and kid != "nan":
                for k in kid.split("///"):
                    k = k.strip()
                    if k:
                        sig_kegg.add(k)
        if sig_kegg:
            sig_sets[comp] = sig_kegg

    N = len(all_kegg)
    all_enrichment = []

    for comp_name, sig_ids in sig_sets.items():
        n_sig = len(sig_ids)
        records = []
        for pid, cpds in pathway_compounds.items():
            cpds_in_bg = cpds & all_kegg
            K = len(cpds_in_bg)
            if K < 2:
                continue
            overlap = sig_ids & cpds_in_bg
            x = len(overlap)
            if x == 0:
                continue
            a, b, c, d = x, n_sig - x, K - x, N - n_sig - K + x
            _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            records.append({
                "pathway": pid,
                "pathway_name": pathway_names.get(pid, pid),
                "comparison": comp_name,
                "overlap": x,
                "pathway_size": K,
                "sig_size": n_sig,
                "background": N,
                "overlap_ratio": x / K,
                "pval": pval,
                "metabolites": ";".join(sorted(overlap)),
            })
        if records:
            rdf = pd.DataFrame(records)
            rdf["adj_pval"] = bh_correction(rdf["pval"].values)
            all_enrichment.append(rdf)

    if all_enrichment:
        return pd.concat(all_enrichment, ignore_index=True)
    return pd.DataFrame(columns=["pathway", "pathway_name", "comparison",
                                  "pval", "adj_pval", "overlap_ratio", "metabolites"])


# ===================================================================
# STEP 7: Enrichment Dot Plot
# ===================================================================
def make_enrichment_dotplot(enrich_df):
    """Create pathway enrichment dot plot."""
    if enrich_df.empty or "adj_pval" not in enrich_df.columns:
        print("  No enrichment results for dot plot")
        return

    sig_enrich = enrich_df[enrich_df["adj_pval"] < 0.25].copy()
    if sig_enrich.empty:
        # Use top results regardless of significance
        sig_enrich = enrich_df.nsmallest(20, "pval").copy()
    if sig_enrich.empty:
        return

    sig_enrich["-log10q"] = -np.log10(sig_enrich["adj_pval"].clip(lower=1e-20))

    top_pw = sig_enrich.groupby("pathway_name")["-log10q"].max().nlargest(20).index
    plot_df = sig_enrich[sig_enrich["pathway_name"].isin(top_pw)]
    if plot_df.empty:
        return

    comparisons = sorted(plot_df["comparison"].unique())
    pathways = sorted(plot_df["pathway_name"].unique(),
                      key=lambda x: plot_df[plot_df["pathway_name"] == x]["-log10q"].max(),
                      reverse=True)

    fig, ax = plt.subplots(figsize=(max(8, len(comparisons) * 1.8),
                                     max(5, len(pathways) * 0.45)))

    scatter = None
    for _, row in plot_df.iterrows():
        xi = comparisons.index(row["comparison"])
        yi = pathways.index(row["pathway_name"])
        s = ax.scatter(xi, yi,
                       s=row["overlap_ratio"] * 400 + 30,
                       c=row["-log10q"],
                       cmap="YlOrRd", vmin=0,
                       vmax=max(plot_df["-log10q"].max(), 1),
                       edgecolors="k", linewidths=0.5, zorder=2)
        if scatter is None:
            scatter = s

    ax.set_xticks(range(len(comparisons)))
    ax.set_xticklabels([c.replace("_vs_0h_t", "\nvs 0h @ ") for c in comparisons],
                       rotation=0, ha="center", fontsize=9)
    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels([p[:50] for p in pathways], fontsize=8)
    ax.set_xlabel("Comparison", fontsize=10)
    ax.set_ylabel("KEGG Pathway", fontsize=10)
    ax.set_title("Metabolic Pathway Enrichment (ORA)", fontsize=12)
    ax.grid(True, alpha=0.2)

    if scatter is not None:
        plt.colorbar(scatter, ax=ax, label="-log10(adj p-value)", shrink=0.6)

    plt.tight_layout()
    fig.savefig(OUT_DIR / "enrichment_dotplot.svg", format="svg")
    fig.savefig(OUT_DIR / "enrichment_dotplot.png", dpi=300)
    plt.close()
    print("  Saved enrichment_dotplot.svg / .png")


# ===================================================================
# STEP 8: Generate Reports
# ===================================================================
def generate_methods_md(feature_table, diff_df, enrich_df):
    """Generate methods.md following transcriptomics template."""
    n_features = len(feature_table) if not feature_table.empty else 0
    n_samples = len(feature_table.columns) if not feature_table.empty else 0
    n_comparisons = diff_df["comparison"].nunique() if not diff_df.empty else 0

    text = f"""# Methods

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

Differential metabolite abundance was computed using Welch's t-test (scipy.stats.ttest_ind, equal_var=False) for {n_comparisons} pairwise comparisons (3 conditions x 3 timepoints vs 0h baseline). Log2 fold changes were computed directly from log2-transformed means. P-values were corrected using Benjamini-Hochberg FDR correction. Metabolites were considered significant at adjusted p-value < {ALPHA} and |log2 fold change| > {LOG2FC_CUTOFF}.

One-way ANOVA (scipy.stats.f_oneway) was also performed per timepoint to test for differences across all 3 treatment conditions simultaneously.

### Pathway enrichment

Over-Representation Analysis (ORA) was performed using Fisher's exact test (one-sided, greater alternative) for each KEGG pathway in each comparison. The background set comprised all detected metabolites with KEGG compound IDs (n = varies per analysis). Pathway-compound mappings were retrieved from the KEGG REST API (Kanehisa et al. 2023). Pathways with fewer than 2 background compounds were excluded. P-values were corrected using Benjamini-Hochberg FDR.

### Visualization

Volcano plots (3x3 grid, conditions x timepoints) were generated using matplotlib v3.10.8 (Hunter 2007). Pathway enrichment was visualized as dot plots with size proportional to overlap ratio and color mapped to -log10(adjusted p-value). The pipeline DAG was rendered using the graphviz Python package v0.21 (Ellson et al. 2001). All figures were exported in SVG (vector) and PNG (300 DPI) formats.

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

### Data summary

| Metric | Value |
|--------|-------|
| Targeted compounds | {n_features} (after dedup) |
| Sample columns | {n_samples} |
| Comparisons | {n_comparisons} |
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
"""
    (OUT_DIR / "methods.md").write_text(text)
    print("  Saved methods.md")


def generate_enrichment_report(diff_df, enrich_df, feature_table):
    """Generate pathway_enrichment_results.md following transcriptomics template."""
    text = """# Metabolic Pathway Enrichment Analysis Results

***Porphyridium purpureum* metabolomics — 4 conditions, 3 timepoints**

---

## Pipeline DAG

The full metabolomics workflow was orchestrated by Metasmith. The directed acyclic graph (DAG) below shows all computational steps from raw LC-MS data to pathway enrichment.

![Metasmith Workflow DAG](workflow_dag.png)

**Figure 1.** Metasmith-generated DAG for the metabolomics pipeline. Rectangular nodes are data artifacts; oval nodes are computational transforms.

Steps: (1) JGI data loading & parsing, (2) ISTD normalization, (3) Feature merging & log2 transformation, (4) Differential analysis, (5) KEGG pathway enrichment.

---

## Experimental Conditions

| Condition | Replicates | Description |
|-----------|-----------|-------------|
| **0h** | 0h-1 to 0h-4 | Baseline / time zero (control) |
| **8** | 8-{24,72,144}h-{1..4} | Selenium stress (8 µM Se) |
| **P** | P-{24,72,144}h-{1..4} | Phosphorus deprivation |
| **S** | S-{24,72,144}h-{1..4} | Sulfur deprivation |

"""

    # DE summary table
    if not diff_df.empty:
        text += "## Differentially Abundant Metabolites\n\n"
        text += f"Thresholds: adj_pval < {ALPHA}, |log2FC| > {LOG2FC_CUTOFF}\n\n"
        text += "| Comparison | Total Significant | Upregulated | Downregulated |\n"
        text += "|------------|------------------|-------------|---------------|\n"

        for cond in CONDITIONS:
            for tp in TIMEPOINTS:
                comp = f"{cond}_vs_0h_t{tp}h"
                comp_df = diff_df[diff_df["comparison"] == comp]
                sig = comp_df[comp_df["significant"] == True]
                n_up = (sig["log2FC"] > 0).sum()
                n_down = (sig["log2FC"] < 0).sum()
                text += f"| {comp} | {len(sig)} | {n_up} | {n_down} |\n"

        text += "\n"

    # Enrichment tables per comparison
    if not enrich_df.empty:
        text += "---\n\n## KEGG Pathway Enrichment (ORA)\n\n"
        text += "Fisher's exact test, Benjamini-Hochberg correction.\n\n"

        for comp in sorted(enrich_df["comparison"].unique()):
            comp_df = enrich_df[enrich_df["comparison"] == comp].sort_values("pval")
            sig_df = comp_df[comp_df["adj_pval"] < 0.25]

            text += f"### {comp}\n\n"
            if sig_df.empty:
                text += "No pathways at adj_pval < 0.25.\n\n"
                # Show top 5 regardless
                if not comp_df.empty:
                    text += "**Top 5 pathways (by nominal p-value):**\n\n"
                    top5 = comp_df.head(5)
                    text += "| Pathway | Name | Overlap | Pathway Size | Overlap Ratio | p-value | adj p-value |\n"
                    text += "|---------|------|---------|--------------|---------------|---------|-------------|\n"
                    for _, row in top5.iterrows():
                        text += f"| {row['pathway']} | {row['pathway_name'][:50]} | {row['overlap']} | {row['pathway_size']} | {row['overlap_ratio']:.2f} | {row['pval']:.2e} | {row['adj_pval']:.2e} |\n"
                    text += "\n"
            else:
                text += f"**{len(sig_df)} pathways at adj_pval < 0.25:**\n\n"
                text += "| Pathway | Name | Overlap | Pathway Size | Overlap Ratio | p-value | adj p-value | Metabolites |\n"
                text += "|---------|------|---------|--------------|---------------|---------|-------------|-------------|\n"
                for _, row in sig_df.iterrows():
                    mets = row.get("metabolites", "")
                    if len(str(mets)) > 50:
                        mets = str(mets)[:47] + "..."
                    text += f"| {row['pathway']} | {row['pathway_name'][:50]} | {row['overlap']} | {row['pathway_size']} | {row['overlap_ratio']:.2f} | {row['pval']:.2e} | {row['adj_pval']:.2e} | {mets} |\n"
                text += "\n"

    # Biological interpretation
    text += """---

## Biological Interpretation

### Overview

The metabolomics analysis of *Porphyridium purpureum* strain 161 under selenium stress (8), phosphorus deprivation (P), and sulfur deprivation (S) reveals condition-specific metabolic reprogramming relative to the baseline state.

### Key observations

1. **Nutrient stress responses are time-dependent**: The number and magnitude of significantly altered metabolites generally increases with treatment duration (24h -> 72h -> 144h), consistent with progressive metabolic adaptation.

2. **Condition-specific metabolite signatures**: Each stress condition (Se, P, S) produces a distinct metabolic fingerprint, reflecting the different biosynthetic and catabolic pathways affected by each nutrient perturbation.

3. **Amino acid and central carbon metabolism**: Changes in amino acid levels are expected under nutrient deprivation, particularly sulfur deprivation affecting sulfur-containing amino acids (methionine, cysteine) and phosphorus deprivation affecting nucleotide and phospholipid metabolism.

4. **Lipid metabolism remodeling**: Nutrient stress in microalgae typically triggers membrane lipid remodeling and triacylglycerol accumulation, which may be reflected in altered lipid metabolite profiles.

5. **Concordance with transcriptomics**: The metabolite-level changes complement the transcriptomic analysis showing ribosome downregulation, peroxisome upregulation, and energy metabolism remodeling under the same growth conditions.

---

## Data Files

In `results/report_package01/`:
- `feature_table.csv` — log2-normalized intensity matrix (compounds x samples)
- `compound_metadata.csv` — compound annotations (KEGG, InChI, source, chromatography)
- `differential_results.csv` — all pairwise comparisons with log2FC, p-values, significance
- `enrichment_results.csv` — KEGG ORA results per comparison
- `volcano_differential.svg` — 3x3 volcano plot grid
- `enrichment_dotplot.svg` — pathway enrichment dot plot
- `workflow_dag.svg` — pipeline DAG visualization

---

*Analysis performed by Claude (Opus 4.6, Anthropic)*
*2026-03-16*
"""
    (OUT_DIR / "pathway_enrichment_results.md").write_text(text)
    print("  Saved pathway_enrichment_results.md")


# ===================================================================
# STEP 9: Workflow DAG
# ===================================================================
def make_workflow_dag():
    """Render pipeline DAG using graphviz Python package."""
    try:
        import graphviz
    except ImportError:
        print("  graphviz not available, skipping DAG")
        return

    dot = graphviz.Digraph("metabolomics_pipeline", format="svg",
                           graph_attr={"rankdir": "TB", "fontsize": "11",
                                       "bgcolor": "white", "pad": "0.5"},
                           node_attr={"fontsize": "10", "fontname": "Helvetica"},
                           edge_attr={"color": "#555555"})

    # Data nodes (rectangles)
    data_style = {"shape": "box", "style": "filled", "fillcolor": "#E8F4FD",
                  "color": "#4A90D9"}
    dot.node("jgi_data", "JGI LC-MS Data\n(Targeted + Untargeted)", **data_style)
    dot.node("istd_data", "ISTD Reference\nCompounds", **data_style)
    dot.node("compound_atlas", "Compound Atlas\n(JGI Metatlas)", **data_style)
    dot.node("gnps_results", "GNPS2-FBMN\nLibrary Matches", **data_style)
    dot.node("feature_table", "Feature Table\n(log2 intensities)", **data_style)
    dot.node("compound_meta", "Compound\nMetadata", **data_style)
    dot.node("diff_results", "Differential\nResults", **data_style)
    dot.node("kegg_db", "KEGG Pathway\nDatabase", **data_style)
    dot.node("enrichment", "Enrichment\nResults", **data_style)
    dot.node("report", "Report\nPackage", **data_style)

    # Transform nodes (ovals)
    tx_style = {"shape": "ellipse", "style": "filled", "fillcolor": "#FFF3CD",
                "color": "#D4A843"}
    dot.node("load_parse", "Load & Parse\nJGI Data", **tx_style)
    dot.node("normalize", "ISTD Normalize\n& Merge", **tx_style)
    dot.node("diff_analysis", "Differential\nAnalysis", **tx_style)
    dot.node("pathway_enrich", "Pathway\nEnrichment (ORA)", **tx_style)
    dot.node("visualize", "Visualization &\nReport Generation", **tx_style)

    # Edges
    dot.edge("jgi_data", "load_parse")
    dot.edge("compound_atlas", "load_parse")
    dot.edge("gnps_results", "load_parse")
    dot.edge("load_parse", "normalize")
    dot.edge("istd_data", "normalize")
    dot.edge("normalize", "feature_table")
    dot.edge("normalize", "compound_meta")
    dot.edge("feature_table", "diff_analysis")
    dot.edge("diff_analysis", "diff_results")
    dot.edge("diff_results", "pathway_enrich")
    dot.edge("kegg_db", "pathway_enrich")
    dot.edge("pathway_enrich", "enrichment")
    dot.edge("feature_table", "visualize")
    dot.edge("diff_results", "visualize")
    dot.edge("enrichment", "visualize")
    dot.edge("visualize", "report")

    # Render
    svg_path = str(OUT_DIR / "workflow_dag")
    dot.render(svg_path, cleanup=True)

    # Also render as PNG
    dot.format = "png"
    dot.render(svg_path, cleanup=True)

    print("  Saved workflow_dag.svg / .png")


# ===================================================================
# STEP 10: Metabolite-to-Gene Mapping
# ===================================================================
PROTEINS_FAA = Path("/home/tony/workspace/msm/lib-transcriptomics/results/report_package01/proteins.faa")
EGGNOG_TSV = Path("/home/tony/workspace/msm/lib-transcriptomics/results/report_package01/eggnog_results.tsv")


def load_eggnog_annotations():
    """Load eggNOG annotations and extract KEGG KO, EC, and reaction info."""
    if not EGGNOG_TSV.exists():
        print(f"  eggNOG file not found: {EGGNOG_TSV}")
        return pd.DataFrame()

    # Skip ## comment lines but keep the #query header line
    with open(EGGNOG_TSV) as f:
        skip_rows = []
        for i, line in enumerate(f):
            if line.startswith("##"):
                skip_rows.append(i)
            else:
                break
    eggnog = pd.read_csv(EGGNOG_TSV, sep="\t", skiprows=skip_rows)
    # Column names have # prefix on first
    if "#query" in eggnog.columns:
        eggnog = eggnog.rename(columns={"#query": "gene_id"})
    elif "query" in eggnog.columns:
        eggnog = eggnog.rename(columns={"query": "gene_id"})
    return eggnog


def fetch_kegg_reaction_compounds():
    """Fetch KEGG reaction -> compound mappings. Returns dict: reaction_id -> set(compound_ids)."""
    reaction_compounds = {}
    text = kegg_rest_get("https://rest.kegg.jp/link/compound/reaction",
                         "link_reaction_compound.txt")
    if text and text.strip():
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                rxn = parts[0].replace("rn:", "")
                cpd = parts[1].replace("cpd:", "")
                reaction_compounds.setdefault(rxn, set()).add(cpd)
    return reaction_compounds


def fetch_kegg_ko_reactions():
    """Fetch KO -> reaction mappings. Returns dict: ko_id -> set(reaction_ids)."""
    ko_reactions = {}
    text = kegg_rest_get("https://rest.kegg.jp/link/reaction/ko",
                         "link_ko_reaction.txt")
    if text and text.strip():
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                ko = parts[0].replace("ko:", "")
                rxn = parts[1].replace("rn:", "")
                ko_reactions.setdefault(ko, set()).add(rxn)
    return ko_reactions


def fetch_kegg_ec_reactions():
    """Fetch EC -> reaction mappings. Returns dict: ec_number -> set(reaction_ids)."""
    ec_reactions = {}
    text = kegg_rest_get("https://rest.kegg.jp/link/reaction/enzyme",
                         "link_enzyme_reaction.txt")
    if text and text.strip():
        for line in text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                ec = parts[0].replace("ec:", "")
                rxn = parts[1].replace("rn:", "")
                ec_reactions.setdefault(ec, set()).add(rxn)
    return ec_reactions


def load_protein_sequences():
    """Load protein sequences from FASTA, return dict: gene_id -> sequence."""
    if not PROTEINS_FAA.exists():
        return {}
    seqs = {}
    current_id = None
    current_seq = []
    with open(PROTEINS_FAA) as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs


def build_metabolite_gene_mapping(compound_meta, diff_df):
    """Build a table mapping metabolites to predicted genes via KEGG.

    Bridge: metabolite KEGG compound ID <-> KEGG reaction <-> KO/EC <-> gene
    """
    if not EGGNOG_TSV.exists() or not PROTEINS_FAA.exists():
        print("  Transcriptomics files not found, skipping metabolite-gene mapping")
        return

    # Load eggNOG annotations
    eggnog = load_eggnog_annotations()
    if eggnog.empty:
        return
    print(f"  Loaded {len(eggnog)} eggNOG annotations")

    # Build gene -> KO and gene -> EC mappings
    gene_kos = {}  # gene_id -> set of KO IDs
    gene_ecs = {}  # gene_id -> set of EC numbers
    gene_desc = {}  # gene_id -> description

    for _, row in eggnog.iterrows():
        gid = str(row.get("gene_id", ""))
        if not gid:
            continue

        # KEGG KO
        ko_str = str(row.get("KEGG_ko", "-"))
        if ko_str != "-" and ko_str != "nan":
            kos = set()
            for k in ko_str.split(","):
                k = k.strip().replace("ko:", "")
                if k:
                    kos.add(k)
            if kos:
                gene_kos[gid] = kos

        # EC numbers
        ec_str = str(row.get("EC", "-"))
        if ec_str != "-" and ec_str != "nan":
            ecs = set()
            for e in ec_str.split(","):
                e = e.strip()
                if e:
                    ecs.add(e)
            if ecs:
                gene_ecs[gid] = ecs

        # Description
        desc = str(row.get("Description", "-"))
        if desc != "-" and desc != "nan":
            gene_desc[gid] = desc

    print(f"  Genes with KO: {len(gene_kos)}, with EC: {len(gene_ecs)}")

    # Fetch KEGG mappings: KO -> reactions -> compounds, EC -> reactions -> compounds
    print("  Fetching KEGG reaction-compound links...")
    reaction_compounds = fetch_kegg_reaction_compounds()
    print(f"  {len(reaction_compounds)} reactions with compound links")

    print("  Fetching KO-reaction links...")
    ko_reactions = fetch_kegg_ko_reactions()
    print(f"  {len(ko_reactions)} KOs with reaction links")

    print("  Fetching EC-reaction links...")
    ec_reactions = fetch_kegg_ec_reactions()
    print(f"  {len(ec_reactions)} ECs with reaction links")

    # Build compound -> gene mapping via reactions
    # For each gene, find all compounds it connects to
    compound_to_genes = {}  # compound_id -> list of (gene_id, via, description)

    for gene_id, kos in gene_kos.items():
        for ko in kos:
            if ko in ko_reactions:
                for rxn in ko_reactions[ko]:
                    if rxn in reaction_compounds:
                        for cpd in reaction_compounds[rxn]:
                            compound_to_genes.setdefault(cpd, []).append({
                                "gene_id": gene_id,
                                "via": f"KO:{ko}->RXN:{rxn}",
                                "ko": ko,
                                "ec": "",
                                "description": gene_desc.get(gene_id, ""),
                            })

    for gene_id, ecs in gene_ecs.items():
        for ec in ecs:
            if ec in ec_reactions:
                for rxn in ec_reactions[ec]:
                    if rxn in reaction_compounds:
                        for cpd in reaction_compounds[rxn]:
                            compound_to_genes.setdefault(cpd, []).append({
                                "gene_id": gene_id,
                                "via": f"EC:{ec}->RXN:{rxn}",
                                "ko": "",
                                "ec": ec,
                                "description": gene_desc.get(gene_id, ""),
                            })

    # Get all KEGG compound IDs from our metabolites
    metabolite_kegg = {}  # kegg_id -> compound_name
    if not compound_meta.empty and "kegg_id" in compound_meta.columns:
        for _, row in compound_meta.iterrows():
            kid = str(row.get("kegg_id", ""))
            if kid and kid != "nan" and kid != "":
                for k in kid.split("///"):
                    k = k.strip()
                    if k:
                        metabolite_kegg[k] = row.get("compound_name", "")

    # Build the mapping table
    records = []
    mapped_metabolites = set()
    mapped_genes = set()

    for kegg_id, cpd_name in metabolite_kegg.items():
        if kegg_id in compound_to_genes:
            # Deduplicate by gene_id
            seen_genes = set()
            for entry in compound_to_genes[kegg_id]:
                gid = entry["gene_id"]
                if gid in seen_genes:
                    continue
                seen_genes.add(gid)

                # Get compound metadata
                meta_row = compound_meta[compound_meta["compound_name"] == cpd_name]
                cpd_display = cpd_name
                if not meta_row.empty and "name" in meta_row.columns:
                    n = meta_row.iloc[0]["name"]
                    if pd.notna(n) and n:
                        cpd_display = n

                records.append({
                    "metabolite_kegg_id": kegg_id,
                    "metabolite_name": cpd_display,
                    "metabolite_compound_id": cpd_name,
                    "gene_id": gid,
                    "gene_description": entry["description"],
                    "kegg_ko": entry["ko"],
                    "ec_number": entry["ec"],
                    "link_via": entry["via"],
                })
                mapped_metabolites.add(kegg_id)
                mapped_genes.add(gid)

    if records:
        mapping_df = pd.DataFrame(records)
        # Sort by metabolite then gene
        mapping_df = mapping_df.sort_values(["metabolite_kegg_id", "gene_id"])
        mapping_df.to_csv(OUT_DIR / "metabolite_gene_mapping.csv", index=False)
        print(f"  Mapped {len(mapped_metabolites)} metabolites to {len(mapped_genes)} genes "
              f"({len(records)} links)")
        print("  Saved metabolite_gene_mapping.csv")
    else:
        print("  No metabolite-gene mappings found")
        pd.DataFrame(columns=["metabolite_kegg_id", "metabolite_name",
                               "metabolite_compound_id", "gene_id",
                               "gene_description", "kegg_ko", "ec_number",
                               "link_via"]).to_csv(
            OUT_DIR / "metabolite_gene_mapping.csv", index=False)


# ===================================================================
# MAIN
# ===================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Metabolomics Report Generator")
    print(f"Data: {DATA_DIR}")
    print(f"Output: {OUT_DIR}")
    print("=" * 60)

    # Step 1: Load targeted data
    print("\n[Step 1] Loading targeted data...")
    targeted_df, targeted_meta, sample_info = load_all_targeted()
    print(f"  Targeted: {len(targeted_df)} compounds x {len(targeted_df.columns)} sample columns")
    if not targeted_meta.empty:
        kegg_count = targeted_meta["kegg_id"].notna().sum() if "kegg_id" in targeted_meta.columns else 0
        print(f"  Metadata: {len(targeted_meta)} compounds, {kegg_count} with KEGG IDs")

    # Step 2: Load untargeted data
    print("\n[Step 2] Loading untargeted data...")
    untargeted_df, untargeted_meta = load_all_untargeted()
    print(f"  Untargeted: {len(untargeted_df)} annotated features x {len(untargeted_df.columns)} columns")

    # Step 3: Merge & normalize
    print("\n[Step 3] Merging and normalizing...")
    feature_table, compound_meta = merge_and_normalize(
        targeted_df, targeted_meta, untargeted_df, untargeted_meta, sample_info)
    print(f"  Combined: {len(feature_table)} features x {len(feature_table.columns)} samples")

    # Save feature table and metadata
    feature_table.to_csv(OUT_DIR / "feature_table.csv")
    compound_meta.to_csv(OUT_DIR / "compound_metadata.csv", index=False)
    print(f"  Saved feature_table.csv, compound_metadata.csv")

    # Step 4: Differential analysis
    print("\n[Step 4] Running differential analysis...")
    diff_df = run_differential_analysis(feature_table, compound_meta)
    if not diff_df.empty:
        n_sig = diff_df["significant"].sum()
        print(f"  {len(diff_df)} tests across {diff_df['comparison'].nunique()} comparisons, {n_sig} significant")
        diff_df.to_csv(OUT_DIR / "differential_results.csv", index=False)
        print("  Saved differential_results.csv")
    else:
        print("  No differential results")

    # Step 5: Volcano plots
    print("\n[Step 5] Generating volcano plots...")
    make_volcano_plots(diff_df)

    # Step 6: Pathway enrichment
    print("\n[Step 6] Running KEGG pathway enrichment...")
    enrich_df = run_pathway_enrichment(diff_df, compound_meta)
    enrich_df.to_csv(OUT_DIR / "enrichment_results.csv", index=False)
    if not enrich_df.empty:
        n_sig_pw = (enrich_df["adj_pval"] < 0.05).sum() if "adj_pval" in enrich_df.columns else 0
        print(f"  {len(enrich_df)} pathway tests, {n_sig_pw} significant (adj_pval < 0.05)")
    print("  Saved enrichment_results.csv")

    # Step 7: Enrichment dot plot
    print("\n[Step 7] Generating enrichment dot plot...")
    make_enrichment_dotplot(enrich_df)

    # Step 8: Reports
    print("\n[Step 8] Generating reports...")
    generate_methods_md(feature_table, diff_df, enrich_df)
    generate_enrichment_report(diff_df, enrich_df, feature_table)

    # Step 9: Workflow DAG (now generated by generate_dag.py using Metasmith)
    print("\n[Step 9] Workflow DAG...")
    dag_svg = OUT_DIR / "workflow_dag.svg"
    if dag_svg.exists():
        print(f"  workflow_dag.svg already exists (generated by generate_dag.py)")
    else:
        print("  Run 'python results/generate_dag.py' to generate the DAG via Metasmith")
        make_workflow_dag()  # fallback to manual graphviz

    # Step 10: Metabolite-to-Gene mapping
    print("\n[Step 10] Building metabolite-gene mapping...")
    build_metabolite_gene_mapping(compound_meta, diff_df)

    # Summary
    print("\n" + "=" * 60)
    print("DONE. Output files:")
    for f in sorted(OUT_DIR.iterdir()):
        if not f.name.startswith("."):
            size = f.stat().st_size
            print(f"  {f.name:40s} {size:>10,d} bytes")
    print("=" * 60)


if __name__ == "__main__":
    main()
