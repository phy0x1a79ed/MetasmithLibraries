"""ECSPr reaction-network pipeline -- evidence weights, host base graphs, and the
observed per-fosmid dR_eff / dI_eff solve.

Ports the scadc `04_reaction_network` observed-perturbation chain into one
importable + CLI module, on top of the verified `ecspr_solver` (SMW/Woodbury):

  * evidence weights          <- 10_evidence_weights.py (_nomination_contributions,
                                 compute_E) -- belief-conservation allocation.
  * per-unit addition weights <- 13_base_graphs.per_unit_weights (per fosmid contig
                                 / per metaG ORF, NOT the source aggregate).
  * host base graphs          <- 13_base_graphs.build_base_graph (EPI300-induced
                                 element-X backbone, edge w_X = atom_count*E_epi300).
  * observed solve            <- 14r_reff_axes.py (per element one SMWGraphContext,
                                 per testable axis a solver reusing the dense inverse,
                                 per fosmid reff_summary) -> reff_axes_report.tsv.
  * ieff derivation           <- reff/06_derive_ieff.py (dI = 1/r_aug - 1/r_base).

Significance (empirical-CCDF vs the reused frozen null) and ablation (LOO +
gene-solo) live alongside in the ablation/significance entry points (added
separately). Reference artifacts consumed as staged inputs: the element bipartite
universe graphs mnx_bipartite_{X}.pkl, the EPI300 base_{X}.pkl (host backbone,
fosmid-independent), and the biomass axes JSONs. Only the fosmid side is computed
fresh from the new inserts' evidence table.

Env: numpy + pandas + scipy + networkx (CPU); torch-CUDA only if --device cuda.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ecspr_solver as _es
from ecspr_solver import (SMWGraphContext, SMWSolver, build_ar2m, reff_base,
                          reff_summary, derive_ieff, DirectSolver, reff_direct,
                          _norm_ar2m)

ELEMENTS = ["C", "N", "S", "P"]
REFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "null",
             "iter", "delta_reff", "r_base", "r_aug", "n_added_rxns"]
IEFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "null",
             "iter", "delta_ieff", "g_base", "g_aug", "n_added_rxns"]


# =====================================================================
# Evidence weights (port of 10_evidence_weights)
# =====================================================================

def nomination_contributions(df: pd.DataFrame) -> pd.DataFrame:
    """Per-row contribution = (w_n / F_n) / L_orf, where the nomination unit is
    (orf, channel, intermediate_id), w_n splits raw_score within (orf, channel),
    F_n = distinct-mnxr fanout, and L_orf = distinct channels for the ORF. Each
    ORF's contributions sum to 1.0 (belief conservation), making leave-one-out an
    exact subtraction downstream."""
    d = df.drop_duplicates(["orf", "channel", "intermediate_id", "mnxr"]).copy()
    nom = (d.groupby(["orf", "channel", "intermediate_id"], sort=False)
           .agg(s_n=("raw_score", "max"), F_n=("mnxr", "nunique"))
           .reset_index())
    grp = nom.groupby(["orf", "channel"], sort=False)
    nom["sum_s"] = grp["s_n"].transform("sum")
    nom["n_nom"] = grp["intermediate_id"].transform("count")
    nonzero = nom["sum_s"] > 0
    nom["w_n"] = np.where(nonzero, nom["s_n"] / nom["sum_s"], 1.0 / nom["n_nom"])
    nom["contrib"] = nom["w_n"] / nom["F_n"]
    out = d.merge(nom[["orf", "channel", "intermediate_id", "contrib"]],
                  on=["orf", "channel", "intermediate_id"], how="left")
    n_lanes = out.groupby("orf")["channel"].transform("nunique")
    out["contrib"] = out["contrib"] / n_lanes
    return out


def _assert_conservation(rows: pd.DataFrame, label: str) -> None:
    per_orf = rows.groupby("orf")["contrib"].sum()
    worst = float((per_orf - 1.0).abs().max()) if len(per_orf) else 0.0
    assert worst < 1e-9, f"{label}: per-orf conservation off by {worst:.2e}"


def compute_E(df_src: pd.DataFrame, label: str = "") -> pd.Series:
    rows = nomination_contributions(df_src)
    _assert_conservation(rows, label or "E")
    return rows.groupby("mnxr")["contrib"].sum()


def compute_weights(ev: pd.DataFrame) -> pd.DataFrame:
    """evidence_weights.parquet: per (source, mnxr) E_full/E_dlec/n_orf."""
    parts = []
    for source, g in ev.groupby("source"):
        e_full = compute_E(g, f"{source}/full")
        dlec = g[g["channel"] == "dl_ec"]
        e_dlec = compute_E(dlec, f"{source}/dlec") if len(dlec) else pd.Series(dtype=float)
        n_orf = g.groupby("mnxr")["orf"].nunique()
        part = (pd.DataFrame({"E_full": e_full})
                .join(e_dlec.rename("E_dlec"), how="outer")
                .join(n_orf.rename("n_orf"), how="left"))
        part["E_full"] = part["E_full"].fillna(0.0)
        part["E_dlec"] = part["E_dlec"].fillna(0.0)
        part["n_orf"] = part["n_orf"].fillna(0).astype(int)
        part.insert(0, "source", source)
        part = part.reset_index().rename(columns={"index": "mnxr"})
        parts.append(part)
    out = pd.concat(parts, ignore_index=True)
    return out[["source", "mnxr", "E_full", "E_dlec", "n_orf"]]


def per_unit_weights(df_src: pd.DataFrame, unit_col: str) -> dict:
    """{unit: {mnxr: E}} with per-ORF-normalized belief, re-aggregated by
    'contig' (per fosmid) or 'orf' (per metaG ORF)."""
    rows = nomination_contributions(df_src)
    if unit_col == "contig":
        rows = rows.copy()
        rows["contig"] = rows["orf"].str.rsplit("_", n=1).str[0]
    agg = rows.groupby([unit_col, "mnxr"], sort=False)["contrib"].sum()
    out: dict = defaultdict(dict)
    for (unit, mnxr), c in agg.items():
        if c > 0:
            out[unit][mnxr] = float(c)
    return dict(out)


# =====================================================================
# Element bipartite + host base graphs (port of _graphs / 13_base_graphs)
# =====================================================================

def load_bipartite(path: Path, element: str) -> nx.Graph:
    """Load the universe element bipartite pkl and drop edges with w_X <= 0."""
    with open(path, "rb") as fh:
        G = pickle.load(fh)
    wk = f"w_{element}"
    drop = [(u, v) for u, v, d in G.edges(data=True) if d.get(wk, 0) <= 0]
    G.remove_edges_from(drop)
    return G


def build_base_graph(element: str, epi_rxns: set, e_epi: dict, G_X: nx.Graph):
    """EPI300-induced element-X base, edge w_X = atom_count * E_epi300(rxn).
    Returns (base_graph, base_met_set)."""
    wk = f"w_{element}"
    base = nx.Graph()
    for r in epi_rxns:
        rn = ("rxn", r)
        if rn not in G_X:
            continue
        er = e_epi.get(r, 0.0)
        if er <= 0:
            continue
        for _, m, d in G_X.edges(rn, data=True):
            atom = d.get(wk, 0.0)
            if atom <= 0:
                continue
            base.add_edge(rn, m, **{wk: float(atom) * float(er)})
    lcc = max(nx.connected_components(base), key=len) if base.number_of_nodes() else set()
    base_mets = {n[1] for n in lcc if n[0] == "met"}
    return base, base_mets


# =====================================================================
# CLI subcommands
# =====================================================================

def cmd_weights(args):
    ev = pd.read_parquet(args.evidence)
    out = compute_weights(ev)
    out.to_parquet(args.out, index=False)
    print(f"[weights] wrote {args.out} ({len(out):,} rows, "
          f"{out['source'].nunique()} sources)")


def cmd_addition_weights(args):
    ev = pd.read_parquet(args.evidence)
    g = ev[ev["source"] == args.source]
    if not len(g):
        raise SystemExit(f"no evidence rows for source={args.source}")
    w = per_unit_weights(g, args.unit)
    with open(args.out, "wb") as fh:
        pickle.dump(w, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[addition-weights] {args.source}/{args.unit}: {len(w):,} units -> {args.out}")


def cmd_base_graphs(args):
    ev = pd.read_parquet(args.evidence)
    epi_rxns = set(ev[ev["source"] == "epi300"]["mnxr"])
    w = pd.read_parquet(args.weights)
    e_epi = {r.mnxr: r.E_full for r in w[w.source == "epi300"].itertuples()}
    axes = json.loads(Path(args.axes).read_text()) if args.axes else {}
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    testable = {}
    for X in args.elements:
        G_X = load_bipartite(Path(args.bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
        base, base_mets = build_base_graph(X, epi_rxns, e_epi, G_X)
        with open(out_dir / f"base_{X}.pkl", "wb") as fh:
            pickle.dump(base, fh, protocol=pickle.HIGHEST_PROTOCOL)
        if axes:
            axX = {k: v for k, v in axes.items() if k.endswith(f"__{X}")}
            testable[X] = [k for k, a in axX.items()
                           if a["source"][0] in base_mets and a["sink"][0] in base_mets]
        print(f"[base-graphs] base_{X}: {base.number_of_nodes():,} nodes / "
              f"{base.number_of_edges():,} edges"
              + (f"; testable {len(testable[X])}/{len(axX)}" if axes else ""))
    if axes:
        (out_dir / "axes_testable.json").write_text(json.dumps(testable, indent=1))


# =====================================================================
# Direct-GEM base graphs (Network A) -- crosswalk + induce, uniform E=1.0
#   Ports scadc validation 59_build_crosswalk.py + 60_build_netA_base.py.
# =====================================================================

def _load_reac_xref(reac_xref: Path):
    """reac_xref.tsv -> (bigg.reaction -> MNXR, kegg.reaction -> MNXR) first-wins."""
    bigg2m, kegg2m = {}, {}
    with open(reac_xref) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1].startswith("MNXR"):
                continue
            if p[0].startswith("bigg.reaction:"):
                bigg2m.setdefault(p[0].split(":", 1)[1], p[1])
            elif p[0].startswith("kegg.reaction:"):
                kegg2m.setdefault(p[0].split(":", 1)[1], p[1])
    return bigg2m, kegg2m


def _universe_sets(bipartite_dir: Path, elements):
    """Per-universe: (MNXR carrying a w_X>0 edge on ANY element, all MNXR rxn nodes)."""
    wedge, nodes = set(), set()
    for X in elements:
        G = load_bipartite(bipartite_dir / f"mnx_bipartite_{X}.pkl", X)
        for n in G.nodes:
            if n[0] != "rxn":
                continue
            nodes.add(n[1])
            if G.degree(n) > 0:          # load_bipartite already dropped w_X<=0 edges
                wedge.add(n[1])
    return wedge, nodes


def gem_crosswalk(model_json: Path, reac_xref: Path, bipartite_dir: Path, elements):
    """GEM reaction -> ONE canonical CURRENT-MNXR, preferring an atom-mapped
    (w_X>0) universe candidate, then source priority bigg > kegg > embedded.
    Returns (DataFrame[rxn_id, mnxr, source, in_universe], stats dict)."""
    import cobra                                   # lazy: only the GEM path needs it
    bigg2m, kegg2m = _load_reac_xref(reac_xref)
    wedge, nodes = _universe_sets(bipartite_dir, elements)

    def as_list(v):
        return [] if not v else (v if isinstance(v, list) else [v])

    m = cobra.io.load_json_model(str(model_json))
    rows, n_wedge, n_node, n_gap, n_unresolved = [], 0, 0, 0, 0
    for r in m.reactions:
        ann = r.annotation or {}
        cands = []                                  # (priority, source, mnxr)
        if r.id in bigg2m:
            cands.append((0, "bigg", bigg2m[r.id]))
        for kk in as_list(ann.get("kegg.reaction")):
            if kk in kegg2m:
                cands.append((1, "kegg", kegg2m[kk]))
        for xx in as_list(ann.get("metanetx.reaction")):
            cands.append((2, "embedded", xx))
        if not cands:
            n_unresolved += 1
            continue
        wc = [c for c in cands if c[2] in wedge]
        nc = [c for c in cands if c[2] in nodes]
        pick = min(wc or nc or cands, key=lambda c: c[0])
        in_uni = pick[2] in wedge
        rows.append(dict(rxn_id=r.id, mnxr=pick[2], source=pick[1], in_universe=in_uni))
        n_wedge += in_uni
        n_node += (not in_uni) and (pick[2] in nodes)
        n_gap += (not in_uni) and (pick[2] not in nodes)
    df = pd.DataFrame(rows).drop_duplicates("rxn_id")
    stats = dict(n_reactions=len(m.reactions), n_resolved=len(df), n_unresolved=n_unresolved,
                 n_wedge=n_wedge, n_node=n_node, n_gap=n_gap)
    return df, stats


def cmd_gem_base_graphs(args):
    """Network A base: induce the atom-mapped universe on the GEM reactome with
    uniform conductance E=1.0 (w_X = pure universe atom-transit count). Reuse-only
    AAM: every GEM MNXR with a w_X>0 universe edge is a base node; MNXR absent from
    the universe are the AAM gap (reported, not mapped here)."""
    elements = args.elements
    bip = Path(args.bipartite_dir)
    xw, stats = gem_crosswalk(Path(args.model), Path(args.reac_xref), bip, elements)
    mnxrs = set(xw.mnxr.dropna().unique())
    e_uniform = {r: 1.0 for r in mnxrs}
    axes = json.loads(Path(args.axes).read_text()) if args.axes else {}
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    xw.to_parquet(out_dir / "gem_rxn_to_mnxr.parquet", index=False)

    testable, union_present = {}, set()
    rep = [f"# Network A (direct GEM) base -- {Path(args.model).stem}\n",
           f"- crosswalk: {stats['n_resolved']}/{stats['n_reactions']} reactions -> "
           f"{len(mnxrs)} unique MNXR (unresolved {stats['n_unresolved']}); uniform E=1.0",
           f"- atom-mapped in universe (base-usable): {stats['n_wedge']}; "
           f"universe-node-but-stripped: {stats['n_node']}; absent (AAM gap): {stats['n_gap']}\n",
           "| element | nodes | edges | LCC met | MNXR reused | axes testable |",
           "|---|---|---|---|---|---|"]
    for X in elements:
        G_X = load_bipartite(bip / f"mnx_bipartite_{X}.pkl", X)
        base, base_mets = build_base_graph(X, mnxrs, e_uniform, G_X)
        present = {n[1] for n in base.nodes if n[0] == "rxn"}
        union_present |= present
        with open(out_dir / f"base_{X}.pkl", "wb") as fh:
            pickle.dump(base, fh, protocol=pickle.HIGHEST_PROTOCOL)
        n_ax = 0
        if axes:
            axX = {k: v for k, v in axes.items() if k.endswith(f"__{X}")}
            testable[X] = [k for k, a in axX.items()
                           if a["source"][0] in base_mets and a["sink"][0] in base_mets]
            n_ax = len(axX)
        rep.append(f"| {X} | {base.number_of_nodes():,} | {base.number_of_edges():,} | "
                   f"{len(base_mets):,} | {len(present):,} | "
                   f"{len(testable.get(X, []))}/{n_ax} |")
        print(f"[gem-base-graphs] base_{X}: {base.number_of_nodes():,} nodes / "
              f"{base.number_of_edges():,} edges; LCC met {len(base_mets):,}; "
              f"MNXR reused {len(present):,}"
              + (f"; testable {len(testable[X])}/{n_ax}" if axes else ""), flush=True)
    if axes:
        (out_dir / "axes_testable.json").write_text(json.dumps(testable, indent=1))
    gap = sorted(mnxrs - union_present)
    rep.append(f"\n**AAM reuse:** {len(union_present)}/{len(mnxrs)} MNXR atom-mapped in the "
               f"universe (reused, no new mapping). **AAM gap:** {len(gap)} MNXR with no w_X>0 "
               f"edge on any element (transport/exchange/polymer or unmapped).\n")
    (out_dir / "netA_coverage.md").write_text("\n".join(rep) + "\n")
    print(f"[gem-base-graphs] AAM reuse {len(union_present)}/{len(mnxrs)}; gap {len(gap)} "
          f"-> {out_dir}")


def _write_tsv(path, cols, rows):
    with open(path, "w", buffering=1) as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


def _solve_element_direct(X, base, ctx, fos_ar2m, axes, ax_ids):
    """Every (fosmid, axis) for one element: one factorization per fosmid.

    Fosmid-OUTER, axis-inner -- the inverse of the Woodbury path's loop, and the
    whole point. `Delta` and `Z_ee` never depended on the axis; only `z` did. The
    old order therefore rebuilt the entire update per axis, which on CPU means p
    sparse solves per (axis, fosmid). Here a fosmid is factored once and every axis
    is a back-solve column against that one LU.

    Results are keyed and emitted axis-outer afterwards, so the output row order is
    unchanged and a diff against the previous table stays readable.
    """
    idx = ctx.idx
    pairs, keep_ax = [], []
    for ax_id in ax_ids:
        ax = axes[ax_id]
        src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
        if src in idx and snk in idx and idx[src] != idx[snk]:
            pairs.append((idx[src], idx[snk]))
            keep_ax.append(ax_id)
    if not pairs:
        return {}, []

    dsolver = DirectSolver(ctx)
    rb_all = dsolver.base_reff(pairs)

    out = {}
    for contig, ar2m in fos_ar2m.items():
        if not ar2m:
            continue
        ar2m_v, _, _ = _norm_ar2m(ctx, ar2m)     # ctx supplies .idx; that is all it reads
        res = reff_direct(dsolver, ar2m_v, pairs, rb_all)
        na = len(ar2m_v)
        for ax_id, (d, rb, raug) in zip(keep_ax, res):
            out[(ax_id, contig)] = (d, rb, raug, na)
    return out, keep_ax


def cmd_solve(args):
    axes = json.loads(Path(args.axes).read_text())
    testable = json.loads(Path(args.testable).read_text())
    with open(args.addition, "rb") as fh:
        fos_w = pickle.load(fh)
    reff_rows, ieff_rows = [], []
    for X in args.elements:
        wk = f"w_{X}"
        with open(Path(args.base_dir) / f"base_{X}.pkl", "rb") as fh:
            base = pickle.load(fh)
        t = time.time()
        ctx = SMWGraphContext(base, device=args.device,
                              dtype=getattr(args, "dtype", "float64"),
                              edge_weight_key=wk)
        base_nodes = set(ctx.nodes)
        G_X = load_bipartite(Path(args.bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
        fos_ar2m = {c: build_ar2m(w, G_X, base_nodes, wk, reinforce=True)
                    for c, w in fos_w.items()}
        n_nonempty = sum(1 for a in fos_ar2m.values() if a)
        ax_ids = testable.get(X, [])
        if args.smoke:
            ax_ids = ax_ids[:3]
        print(f"[solve] [{X}] n_LCC={ctx.n} ctx {time.time()-t:.1f}s  "
              f"fosmids-with-additions {n_nonempty}/{len(fos_ar2m)}  axes {len(ax_ids)}  "
              f"path={args.path}", flush=True)

        t = time.time()
        if args.path == "direct":
            got, keep_ax = _solve_element_direct(X, base, ctx, fos_ar2m, axes, ax_ids)
        else:
            got, keep_ax = {}, []
            for ax_id in ax_ids:
                ax = axes[ax_id]
                src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
                try:
                    solver = SMWSolver(base, [src], [snk], edge_weight_key=wk, context=ctx)
                except ValueError:
                    continue
                keep_ax.append(ax_id)
                rb = reff_base(solver)
                for contig, ar2m in fos_ar2m.items():
                    if not ar2m:
                        continue
                    d, _rb, raug, na = reff_summary(solver, ar2m, r_base=rb)
                    got[(ax_id, contig)] = (d, rb, raug, na)
        print(f"[solve] [{X}] {len(got)} cells in {time.time()-t:.1f}s", flush=True)

        for ax_id in keep_ax:
            ax = axes[ax_id]
            for contig in fos_ar2m:
                if (ax_id, contig) not in got:
                    continue
                d, rb, raug, na = got[(ax_id, contig)]
                di, gb, ga = derive_ieff(rb, raug)
                common = dict(element=X, axis_id=ax_id, source_hub=ax["source"][0],
                              sink_hub=ax["sink"][0], fosmid=contig, null="obs",
                              iter=0, n_added_rxns=na)
                reff_rows.append({**common, "delta_reff": d, "r_base": rb, "r_aug": raug})
                ieff_rows.append({**common, "delta_ieff": di, "g_base": gb, "g_aug": ga})
    _write_tsv(args.out_reff, REFF_COLS, reff_rows)
    _write_tsv(args.out_ieff, IEFF_COLS, ieff_rows)
    print(f"[solve] wrote {args.out_reff} + {args.out_ieff} ({len(reff_rows)} rows each)")


def cmd_derive_ieff(args):
    df = pd.read_csv(args.reff, sep="\t")
    di, gb, ga = zip(*[derive_ieff(rb, ra) for rb, ra in zip(df["r_base"], df["r_aug"])]) \
        if len(df) else ([], [], [])
    df["delta_ieff"] = list(di); df["g_base"] = list(gb); df["g_aug"] = list(ga)
    keep = [c for c in IEFF_COLS if c in df.columns]
    df[keep].to_csv(args.out, sep="\t", index=False)
    print(f"[derive-ieff] wrote {args.out} ({len(df)} rows)")


def parse_args():
    ap = argparse.ArgumentParser(description="ECSPr reaction-network pipeline")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("weights"); p.set_defaults(fn=cmd_weights)
    p.add_argument("--evidence", required=True); p.add_argument("--out", required=True)

    p = sub.add_parser("addition-weights"); p.set_defaults(fn=cmd_addition_weights)
    p.add_argument("--evidence", required=True); p.add_argument("--out", required=True)
    p.add_argument("--source", default="fosmid"); p.add_argument("--unit", default="contig")

    p = sub.add_parser("base-graphs"); p.set_defaults(fn=cmd_base_graphs)
    p.add_argument("--evidence", required=True); p.add_argument("--weights", required=True)
    p.add_argument("--bipartite-dir", required=True); p.add_argument("--axes", default=None)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--elements", nargs="+", default=ELEMENTS)

    p = sub.add_parser("gem-base-graphs"); p.set_defaults(fn=cmd_gem_base_graphs)
    p.add_argument("--model", required=True, help="curated GEM JSON (iML1515 / iECDH10B)")
    p.add_argument("--reac-xref", required=True, help="MetaNetX reac_xref.tsv (current)")
    p.add_argument("--bipartite-dir", required=True); p.add_argument("--axes", default=None)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--elements", nargs="+", default=ELEMENTS)

    p = sub.add_parser("solve"); p.set_defaults(fn=cmd_solve)
    p.add_argument("--addition", required=True); p.add_argument("--base-dir", required=True)
    p.add_argument("--bipartite-dir", required=True); p.add_argument("--axes", required=True)
    p.add_argument("--testable", required=True)
    p.add_argument("--out-reff", required=True); p.add_argument("--out-ieff", required=True)
    p.add_argument("--elements", nargs="+", default=ELEMENTS)
    p.add_argument("--device", default="cpu"); p.add_argument("--smoke", action="store_true")
    p.add_argument("--dtype", default="float64", choices=["float64", "float32"],
                   help="solve precision; float64 keeps parity and runs on H100 GPU")
    p.add_argument("--path", default="direct", choices=["direct", "woodbury"],
                   help="direct: one factorization per fosmid, all axes as back-solves "
                        "(~65x on C). woodbury: the incumbent update; kept because it "
                        "is independently verified against a dense rebuild and so "
                        "serves as a referent for direct, sharing none of its math.")

    p = sub.add_parser("derive-ieff"); p.set_defaults(fn=cmd_derive_ieff)
    p.add_argument("--reff", required=True); p.add_argument("--out", required=True)

    return ap.parse_args()


if __name__ == "__main__":
    a = parse_args()
    a.fn(a)
