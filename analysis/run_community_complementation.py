"""Community ECSPr complementation -- which biomass axes each member cannot complete
alone but CAN through the community (distributed metabolism / auxotrophy rescue).

For each organism O and each biomass axis (source precursor -> biomass sink):
  * reach_solo : is snk reachable from src in O's SOLO graph -- O's cytoplasm + O's own
                 transporters to an ENV that contains ONLY O (no partner supply)?
  * reach_full : is it reachable in the full 3-member community graph, where O can now
                 import intermediates its partners export to the shared ENV?
A "community-complemented" axis is reach_full & ~reach_solo: O outsources part of that
biosynthetic route to the community. We then attribute it -- for each partner Pt, does
the pair {O, Pt} alone restore reachability? -- naming WHO supplies the missing step.

Unlike raw conductance (which any transporter-bearing bacterium inflates), a
complementation is a discrete dependence: O structurally needs a partner to reach the
sink. Reported per member and per (dependent -> supplier) direction, foregrounding
Erythrobacter, to substantiate its role in the community.

Env: numpy + pandas + networkx (CPU). Reuses ecspr_community.
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "resources" / "lib"))
sys.path.insert(0, str(HERE))
import ecspr_community as ec  # noqa: E402
from run_community_ecspr import load_axes  # noqa: E402
from ecspr_solver import SMWGraphContext  # noqa: E402

ELEMENTS = ["C", "N", "S", "P"]


def reachable_in(G, ctx_cache, src, snk):
    """Bool: is snk reachable from src in G (same connected component, finite R_eff)?
    Uses a cached connectivity structure keyed by the graph object id."""
    if src not in G or snk not in G or src == snk:
        return False
    cc = ec.nx.node_connected_component(G, src)
    return snk in cc


def run(real_w, real_t, bipartite_dir, axes, out_tsv):
    orgs = sorted(real_w)
    rows = []
    for X in ELEMENTS:
        axX = [a for a in axes if a["element"] == X]
        if not axX:
            continue
        G_X = ec.load_bipartite(Path(bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
        atoms = ec.metabolite_atoms(G_X, X)

        def build(members):
            return ec.build_community_graph(
                X, {o: real_w[o] for o in members},
                {o: real_t.get(o, {}) for o in members}, G_X, atoms=atoms)

        G_full = build(orgs)
        G_solo = {o: build([o]) for o in orgs}
        G_pair = {(o, pt): build([o, pt]) for o in orgs for pt in orgs if o != pt}

        for a in axX:
            for O in orgs:
                src, snk = ("met", a["source"], O), ("met", a["sink"], O)
                r_solo = reachable_in(G_solo[O], None, src, snk)
                r_full = reachable_in(G_full, None, src, snk)
                complemented = bool(r_full and not r_solo)
                suppliers = []
                if complemented:
                    for pt in orgs:
                        if pt == O:
                            continue
                        if reachable_in(G_pair[(O, pt)], None, src, snk):
                            suppliers.append(pt)
                rows.append(dict(
                    element=X, axis_id=a["axis_id"], category=a["category"],
                    source=a["source"], sink=a["sink"], member=O,
                    reach_solo=int(r_solo), reach_full=int(r_full),
                    complemented=int(complemented),
                    suppliers=";".join(suppliers)))
        print(f"[compl] [{X}] {len(axX)} axes x {len(orgs)} members done", flush=True)

    df = pd.DataFrame(rows)
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[compl] wrote {out_tsv} ({len(df):,} rows)")
    _summary(df)
    return df


def _summary(df):
    print("\n[compl] per-member biosynthetic self-sufficiency across biomass axes:")
    for O, g in df.groupby("member"):
        tot = len(g)
        solo = int(g["reach_solo"].sum())
        comp = int(g["complemented"].sum())
        print(f"    {O}: reachable-alone {solo}/{tot}, "
              f"community-complemented {comp} (needs a partner)")
    comp = df[df["complemented"] == 1]
    if comp.empty:
        print("\n[compl] no community-complemented axes (every member self-sufficient "
              "on testable axes).")
        return
    print(f"\n[compl] {len(comp)} community-complemented axes "
          "(member cannot complete alone, community restores):")
    # dependent -> supplier edges
    edges = {}
    for r in comp.itertuples(index=False):
        for s in (r.suppliers.split(";") if r.suppliers else []):
            edges[(r.member, s)] = edges.get((r.member, s), 0) + 1
    print("\n[compl] dependence edges (dependent needs <- supplier), n axes:")
    for (dep, sup), n in sorted(edges.items(), key=lambda x: -x[1]):
        star = "  <-- Erythrobacter" if "ERY" in (dep, sup) else ""
        print(f"    {dep} needs <- {sup}: {n} axes{star}")
    print("\n[compl] complemented axes detail (first 20):")
    for r in comp.head(20).itertuples(index=False):
        print(f"    [{r.element}] {r.member} {r.category:10s} {r.source}->{r.sink}  "
              f"supplied by: {r.suppliers or '(none singly)'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real-weights", required=True)
    ap.add_argument("--real-transporters", required=True)
    ap.add_argument("--bipartite-dir", required=True)
    ap.add_argument("--axes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run(pickle.load(open(a.real_weights, "rb")),
        pickle.load(open(a.real_transporters, "rb")),
        Path(a.bipartite_dir), load_axes(Path(a.axes)), Path(a.out))


if __name__ == "__main__":
    main()
