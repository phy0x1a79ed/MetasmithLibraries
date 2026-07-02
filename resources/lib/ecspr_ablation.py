"""ECSPr per-gene ablation -- leave-one-ORF-out (LOO) + gene-solo, both metrics.

Ports scadc `04_reaction_network/reff/ablation/{_ablation.py,30_reff_ablation.py}`
onto the verified `ecspr_solver` + `ecspr_network`, and ADDS the net-new gene-solo
variant (which did not exist in scadc). For each focal fosmid, on each testable
biomass axis, three readouts per ORF:

  * full insert           dR_eff(full) / dI_eff(full)  (== reff_axes_report.tsv)
  * LOO ("minus a gene")  importance(orf)   = dR_eff(full) - dR_eff(full minus orf)
                          importance_ieff   = 1/r_aug_full - 1/r_aug_loo
  * gene-solo ("gene alone") dR_eff_solo(orf) = dR_eff of just that ORF's reactions
                          injected onto the host; dI_eff_solo likewise.

Exactness: `nomination_contributions` computes each ORF's per-reaction contribution
strictly within that ORF, so E_fosmid[mnxr] = sum_orf contrib(orf, mnxr). LOO is an
exact subtraction (no re-normalization of survivors) and gene-solo is just that
ORF's own map. All three are >= 0 by Rayleigh monotonicity.

Batching: per (element, contig) the [full, loo_0.., solo_0..] addition maps are
built once and reused across the element's axes; each (contig, axis) is one
reff_batch call. Reuses base_{X}.pkl + mnx_bipartite_{X}.pkl + axes JSONs as
staged reference inputs; the fosmid evidence + addition weights are fresh.

Env: numpy + pandas + networkx (CPU); torch-CUDA only if --device cuda.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ecspr_network as en
from ecspr_solver import SMWGraphContext, SMWSolver, build_ar2m, reff_base, reff_summary, REFF_EPS

ELEMENTS = ["C", "N", "S", "P"]
EPS = 1e-12          # r-floor guard, matches derive_ieff
W_EPS = 1e-12        # drop sub-noise per-mnxr weights after subtraction
COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "orf", "r_base",
        "delta_reff_full", "delta_reff_loo", "importance",
        "delta_ieff_full", "delta_ieff_loo", "importance_ieff",
        "delta_reff_solo", "delta_ieff_solo", "n_rxn_orf"]
_CHANNEL_RANK = {"kofam": 0, "dl_ec": 1, "uniref50_dr": 2, "pbert_transfer": 3}


def _g(r: float) -> float:
    return 1.0 / max(r, EPS)


# --- per-ORF decomposition (port of _ablation) ------------------------------

def _fosmid_rows(ev: pd.DataFrame, contigs=None) -> pd.DataFrame:
    fos = ev[ev.source == "fosmid"].copy()
    fos["contig"] = fos["orf"].str.rsplit("_", n=1).str[0]
    if contigs is not None:
        fos = fos[fos.contig.isin(set(contigs))]
    return fos


def per_orf_contribs(ev: pd.DataFrame, contigs=None) -> dict:
    """{contig: {orf: {mnxr: contrib}}} -- within-ORF normalization, so it sums
    (per contig, mnxr) to the fosmid addition map."""
    rows = en.nomination_contributions(_fosmid_rows(ev, contigs))
    rows = rows.copy()
    rows["contig"] = rows["orf"].str.rsplit("_", n=1).str[0]
    agg = rows.groupby(["contig", "orf", "mnxr"], sort=False)["contrib"].sum()
    out: dict = defaultdict(lambda: defaultdict(dict))
    for (contig, orf, mnxr), c in agg.items():
        if c > 0:
            out[contig][orf][mnxr] = float(c)
    return {c: {o: dict(m) for o, m in orfs.items()} for c, orfs in out.items()}


def loo_map(full_map: dict, orf_contrib: dict) -> dict:
    m = dict(full_map)
    for mnxr, c in orf_contrib.items():
        if mnxr in m:
            m[mnxr] -= c
            if m[mnxr] <= W_EPS:
                del m[mnxr]
    return m


def orf_functions(ev: pd.DataFrame, contigs=None) -> dict:
    """{orf: {intermediate_id, intermediate_name, channel, mnxr, raw_score}} --
    most interpretable nomination per ORF (informative channel first)."""
    fos = _fosmid_rows(ev, contigs).copy()
    fos["_crank"] = fos["channel"].map(_CHANNEL_RANK).fillna(9).astype(int)
    fos = fos.sort_values(["_crank", "raw_score"], ascending=[True, False]).drop_duplicates("orf")
    cols = ["intermediate_id", "intermediate_name", "channel", "mnxr", "raw_score"]
    return {r.orf: {c: getattr(r, c) for c in cols} for r in fos.itertuples()}


def verify_full_maps(per_orf: dict, fos_w: dict, tol: float = 1e-9) -> float:
    worst = 0.0
    for contig, orfs in per_orf.items():
        recon: dict = defaultdict(float)
        for m in orfs.values():
            for mnxr, c in m.items():
                recon[mnxr] += c
        shipped = fos_w.get(contig, {})
        for k in set(recon) | set(shipped):
            worst = max(worst, abs(recon.get(k, 0.0) - shipped.get(k, 0.0)))
    assert worst < tol, f"per-ORF sum vs shipped map off by {worst:.2e}"
    return worst


def _reff_batch(solver, maps, batch_size, pair=0):
    """CPU batch = per-map reff_summary loop (GPU path handled by solver)."""
    rb = reff_base(solver, pair)
    return np.array([reff_summary(solver, a, pair=pair, r_base=rb)[0] for a in maps])


def axis_endpoints(ax_id: str, axes: dict):
    a = axes.get(ax_id)
    if a is not None:
        return a["source"][0], a["sink"][0]
    parts = ax_id.split("__")
    if len(parts) != 4:
        raise KeyError(f"cannot resolve endpoints for axis {ax_id!r}")
    return parts[1], parts[2]


def build_maps_for_element(fosmids, fos_w, per_orf, G_X, base_nodes, wk):
    """{contig: (orfs, maps)} with maps = [full] + loo(per orf) + solo(per orf)."""
    out = {}
    for contig in fosmids:
        full_map = fos_w.get(contig, {})
        ar2m_full = build_ar2m(full_map, G_X, base_nodes, wk, reinforce=True)
        if not ar2m_full:
            continue
        orfs = sorted(per_orf.get(contig, {}))
        loo = [build_ar2m(loo_map(full_map, per_orf[contig][o]), G_X, base_nodes, wk, reinforce=True)
               for o in orfs]
        solo = [build_ar2m(per_orf[contig][o], G_X, base_nodes, wk, reinforce=True)
                for o in orfs]
        out[contig] = (orfs, [ar2m_full] + loo + solo)
    return out


def run(args):
    ev = pd.read_parquet(args.evidence)
    with open(args.addition, "rb") as fh:
        fos_w = pickle.load(fh)
    axes = json.loads(Path(args.axes).read_text())
    testable = json.loads(Path(args.testable).read_text())
    fosmids = args.fosmids if args.fosmids else sorted(c for c, m in fos_w.items() if m)
    per_orf = per_orf_contribs(ev, fosmids)
    worst = verify_full_maps(per_orf, fos_w)
    print(f"[ablation] {len(fosmids)} fosmids; per-ORF sum vs shipped: max|d|={worst:.2e}", flush=True)

    fh = open(args.out, "w", buffering=1)
    fh.write("\t".join(COLS) + "\n")
    for X in args.elements:
        wk = f"w_{X}"
        with open(Path(args.base_dir) / f"base_{X}.pkl", "rb") as f:
            base = pickle.load(f)
        t = time.time()
        ctx = SMWGraphContext(base, device=args.device,
                              dtype=("float32" if args.device != "cpu" else "float64"),
                              edge_weight_key=wk)
        base_nodes = set(ctx.nodes)
        G_X = en.load_bipartite(Path(args.bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
        maps_by_contig = build_maps_for_element(fosmids, fos_w, per_orf, G_X, base_nodes, wk)
        ax_ids = testable.get(X, [])
        if args.smoke:
            ax_ids = ax_ids[:3]
        print(f"[ablation] [{X}] n_LCC={ctx.n} ctx {time.time()-t:.1f}s  "
              f"contigs {len(maps_by_contig)}/{len(fosmids)}  axes {len(ax_ids)}", flush=True)
        for ax_id in ax_ids:
            src_m, snk_m = axis_endpoints(ax_id, axes)
            try:
                solver = SMWSolver(base, [("met", src_m)], [("met", snk_m)],
                                   edge_weight_key=wk, context=ctx)
            except ValueError:
                continue
            rb = reff_base(solver)
            g_base = _g(rb)
            for contig, (orfs, maps) in maps_by_contig.items():
                deltas = _reff_batch(solver, maps, args.batch_size)
                n = len(orfs)
                d_full = float(deltas[0])
                g_full = _g(rb - d_full)
                di_full = max(g_full - g_base, 0.0)
                for i, o in enumerate(orfs):
                    d_loo = float(deltas[1 + i])
                    d_solo = float(deltas[1 + n + i])
                    imp = max(d_full - d_loo, 0.0)
                    if imp < REFF_EPS:
                        imp = 0.0
                    g_loo = _g(rb - d_loo)
                    di_loo = max(g_loo - g_base, 0.0)
                    imp_i = max(g_full - g_loo, 0.0)
                    g_solo = _g(rb - d_solo)
                    di_solo = max(g_solo - g_base, 0.0)
                    row = dict(element=X, axis_id=ax_id, source_hub=src_m, sink_hub=snk_m,
                               fosmid=contig, orf=o, r_base=f"{rb:.10g}",
                               delta_reff_full=f"{d_full:.10g}", delta_reff_loo=f"{d_loo:.10g}",
                               importance=f"{imp:.10g}",
                               delta_ieff_full=f"{di_full:.10g}", delta_ieff_loo=f"{di_loo:.10g}",
                               importance_ieff=f"{imp_i:.10g}",
                               delta_reff_solo=f"{d_solo:.10g}", delta_ieff_solo=f"{di_solo:.10g}",
                               n_rxn_orf=len(per_orf[contig][o]))
                    fh.write("\t".join(str(row.get(c, "")) for c in COLS) + "\n")
    fh.close()
    print(f"[ablation] wrote {args.out}", flush=True)


def parse_args():
    ap = argparse.ArgumentParser(description="ECSPr per-gene ablation (LOO + gene-solo)")
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--addition", required=True)
    ap.add_argument("--base-dir", required=True)
    ap.add_argument("--bipartite-dir", required=True)
    ap.add_argument("--axes", required=True)
    ap.add_argument("--testable", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--elements", nargs="+", default=ELEMENTS)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--fosmids", nargs="+", default=None)
    ap.add_argument("--smoke", action="store_true")
    return ap.parse_args()


if __name__ == "__main__":
    run(parse_args())
