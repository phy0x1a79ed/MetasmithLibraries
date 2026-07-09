"""Community ECSPr driver -- detect distributed metabolism across the 3-member
Nostoc/Erythrobacter/Allorhizobium community.

Pipeline:
  1. reaction weights  : per-organism {MNXR: E_r} from the DL functional-annotation
                         evidence table (belief-conserving, ecspr_network).
  2. transporter belief: per-organism {MNXM: belief} (ecspr_transporters, TCDB).
  3. community graph   : per element X, compartmented + transporter-gated
                         (ecspr_community.build_community_graph).
  4. detection         : for each biomass axis (src->snk) and each ORDERED organism
                         pair (A=producer, B=consumer), the cross-organism effective
                         conductance g_comm( src@A -> snk@B ) on the community graph.
                         The disjoint baseline (transport removed) is 0 by construction
                         (proven in the selftest), so any finite g_comm IS a
                         partnership-enabled distributed-metabolism route. We also
                         report the same-organism control g( src@A -> snk@A ) and the
                         full-community g.

Output: community_handoffs.tsv  (element, axis_id, category, source, sink,
        producer, consumer, r_comm, g_comm, testable) + a ranked summary.

Env: numpy + pandas + networkx + scipy (CPU). Reuses ecspr_solver/community/network.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "resources" / "lib"))
import ecspr_community as ec  # noqa: E402
from ecspr_solver import SMWGraphContext, SMWSolver, reff_base  # noqa: E402
from ecspr_network import nomination_contributions  # noqa: E402

ELEMENTS = ["C", "N", "S", "P"]


# =====================================================================
# inputs
# =====================================================================

def org_weights_from_evidence(evidence_parquet: Path) -> dict:
    """{organism: {MNXR: E}} from a DL evidence table carrying an `organism` column,
    using the belief-conserving nomination_contributions allocation (each ORF's
    belief sums to 1.0, split across channels/nominations/MNXR fanout)."""
    ev = pd.read_parquet(evidence_parquet)
    assert "organism" in ev.columns, "evidence table needs an `organism` column"
    # nomination_contributions preserves the input columns (incl. organism) and adds
    # `contrib`; each ORF belongs to exactly one organism so grouping is unambiguous.
    rows = nomination_contributions(ev)
    agg = rows.groupby(["organism", "mnxr"])["contrib"].sum()
    out: dict = {}
    for (org, mnxr), c in agg.items():
        if c > 0:
            out.setdefault(org, {})[mnxr] = float(c)
    return out


def load_axes(axes_testable: Path) -> list:
    """Return [{axis_id, category, source, sink, element}] parsed from the axis ids."""
    data = json.loads(Path(axes_testable).read_text())
    axes = []
    for X, ids in data.items():
        for aid in ids:
            p = ec.parse_axis_id(aid)
            if p:
                axes.append({"axis_id": aid, **p})
    return axes


# =====================================================================
# solve helper -- one shared context per graph, reused across axes
# =====================================================================

def solve_pairs_on_graph(G, element: str, endpoint_pairs: list) -> dict:
    """Given a community graph and a list of (key, src_node, snk_node), return
    {key: R_eff} reusing ONE SMW context (sparse LU) for all pairs. R_eff is +inf if
    an endpoint is absent or outside the graph's largest connected component."""
    wk = f"w_{element}"
    out = {}
    if G.number_of_edges() == 0:
        return {k: float("inf") for k, _, _ in endpoint_pairs}
    ctx = SMWGraphContext(G, device="cpu", edge_weight_key=wk)
    lcc = ctx.lcc
    for key, src, snk in endpoint_pairs:
        if src not in lcc or snk not in lcc or src == snk:
            out[key] = float("inf")
            continue
        try:
            solver = SMWSolver(G, [src], [snk], edge_weight_key=wk, context=ctx)
            out[key] = float(reff_base(solver)) if solver.pairs else float("inf")
        except ValueError:
            out[key] = float("inf")
    return out


# =====================================================================
# detection
# =====================================================================

def detect(org_weights: dict, transporters: dict, bipartite_dir: Path,
           axes: list, out_tsv: Path, currency_set: set | None = None):
    orgs = sorted(org_weights)
    print(f"[detect] organisms: {orgs}")
    ordered_pairs = [(a, b) for a in orgs for b in orgs if a != b]

    rows = []
    for X in ELEMENTS:
        axX = [a for a in axes if a["element"] == X]
        if not axX:
            continue
        G_X = ec.load_bipartite(Path(bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
        atoms = ec.metabolite_atoms(G_X, X)

        # full community graph (all organisms) -- context reused across axes
        G_full = ec.build_community_graph(X, org_weights, transporters, G_X,
                                          currency_set=currency_set, atoms=atoms)
        stats = ec.assert_compartment_isolation(G_full)
        print(f"[detect] [{X}] {len(axX)} axes | full graph "
              f"{stats['n_nodes']:,} nodes / {stats['n_edges']:,} edges / "
              f"{stats['n_transport_edges']:,} transport", flush=True)

        # per-pair community graphs (A+B) and per-organism solo graphs
        pair_graphs = {}
        for (a, b) in ordered_pairs:
            key = frozenset((a, b))
            if key not in pair_graphs:
                pair_graphs[key] = ec.build_community_graph(
                    X, {a: org_weights[a], b: org_weights[b]},
                    {a: transporters.get(a, {}), b: transporters.get(b, {})},
                    G_X, currency_set=currency_set, atoms=atoms)
        solo_graphs = {o: ec.build_community_graph(
            X, {o: org_weights[o]}, {o: transporters.get(o, {})},
            G_X, currency_set=currency_set, atoms=atoms) for o in orgs}

        # assemble endpoint pairs per graph, solve with shared contexts
        # full-community cross-org conductance + same-org control:
        full_pairs = []
        for a in axX:
            for (A, B) in ordered_pairs:
                full_pairs.append(((a["axis_id"], A, B),
                                   ("met", a["source"], A), ("met", a["sink"], B)))
            for O in orgs:  # same-organism control
                full_pairs.append(((a["axis_id"], O, O),
                                   ("met", a["source"], O), ("met", a["sink"], O)))
        full_reff = solve_pairs_on_graph(G_full, X, full_pairs)

        # pairwise-isolated community (A+B only) cross-org conductance
        pair_reff = {}
        for key, Gpair in pair_graphs.items():
            a, b = tuple(sorted(key))
            eps = []
            for ax in axX:
                for (A, B) in [(a, b), (b, a)]:
                    eps.append(((ax["axis_id"], A, B),
                                ("met", ax["source"], A), ("met", ax["sink"], B)))
            pair_reff.update(solve_pairs_on_graph(Gpair, X, eps))

        for a in axX:
            for (A, B) in ordered_pairs:
                r_full = full_reff.get((a["axis_id"], A, B), float("inf"))
                r_pair = pair_reff.get((a["axis_id"], A, B), float("inf"))
                r_ctrlB = full_reff.get((a["axis_id"], B, B), float("inf"))
                g_full = 0.0 if not np.isfinite(r_full) or r_full <= 0 else 1.0 / r_full
                g_pair = 0.0 if not np.isfinite(r_pair) or r_pair <= 0 else 1.0 / r_pair
                rows.append(dict(
                    element=X, axis_id=a["axis_id"], category=a["category"],
                    source=a["source"], sink=a["sink"], producer=A, consumer=B,
                    r_full=r_full, g_full=g_full,
                    r_pair=r_pair, g_pair=g_pair,
                    consumer_solo_reachable=int(np.isfinite(r_ctrlB)),
                ))

    df = pd.DataFrame(rows)
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[detect] wrote {out_tsv} ({len(df):,} rows)")
    _summary(df)
    return df


def _summary(df: pd.DataFrame):
    if df.empty:
        print("[detect] no rows"); return
    fin = df[df["g_pair"] > 0]
    print(f"\n[detect] cross-organism handoffs with finite conductance (pairwise A+B):"
          f" {len(fin):,} / {len(df):,} (axis x ordered-pair)")
    print("\n  by ordered pair (producer -> consumer): n_axes with a route, sum g_pair")
    g = (fin.groupby(["producer", "consumer"])
             .agg(n_axes=("axis_id", "nunique"), sum_g=("g_pair", "sum"))
             .reset_index().sort_values("sum_g", ascending=False))
    for r in g.itertuples(index=False):
        star = "  <-- Erythrobacter" if ("ERY" in (r.producer, r.consumer)) else ""
        print(f"    {r.producer} -> {r.consumer}: {r.n_axes} axes, "
              f"sum_g={r.sum_g:.3f}{star}")
    print("\n  top 15 individual handoffs (pairwise):")
    top = fin.sort_values("g_pair", ascending=False).head(15)
    for r in top.itertuples(index=False):
        print(f"    [{r.element}] {r.producer}->{r.consumer} {r.category:10s} "
              f"{r.source}->{r.sink}  g={r.g_pair:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evidence", required=True, help="DL evidence parquet (organism col)")
    ap.add_argument("--transporters", required=True, help="transporter_belief.pkl {org:{mnxm:belief}}")
    ap.add_argument("--bipartite-dir", required=True)
    ap.add_argument("--axes", required=True, help="axes_testable.json (ids encode src/sink)")
    ap.add_argument("--out", required=True, help="community_handoffs.tsv")
    ap.add_argument("--org-weights-out", default=None, help="optional: dump {org:{mnxr:E}} pickle")
    a = ap.parse_args()

    org_weights = org_weights_from_evidence(Path(a.evidence))
    for o, w in org_weights.items():
        print(f"[detect] {o}: {len(w):,} reactions (sum E {sum(w.values()):.1f})")
    if a.org_weights_out:
        pickle.dump(org_weights, open(a.org_weights_out, "wb"))
    transporters = pickle.load(open(a.transporters, "rb"))
    axes = load_axes(Path(a.axes))
    detect(org_weights, transporters, Path(a.bipartite_dir), axes, Path(a.out))


if __name__ == "__main__":
    main()
