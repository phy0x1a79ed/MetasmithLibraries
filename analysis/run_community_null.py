"""Community ECSPr null model -- is a handoff specific to these community members?

For each biomass axis and each ORDERED real pair (A -> B), the observed cross-organism
conductance g_obs is compared against a "swap-a-member" null: replace A (producer swap)
or B (consumer swap) with each of N random, unrelated MAGs annotated by the IDENTICAL
lanes, and recompute the conductance. The pooled null (2N values) gives an empirical
CCDF p-value

    p = (1 + #{ null >= g_obs }) / (1 + n)          (floored at 1/(n+1))

BH-corrected across all axis x ordered-pair tests. A small q means the real members
route that axis's handoff better than random genomes would -- i.e. the distributed-
metabolism link is specific to the community, not a generic consequence of "any two
bacteria share transporters." This is exactly the manuscript's Erythrobacter question:
does ERY contribute more than a random genome swapped into its place?

Env: numpy + pandas + networkx + scipy (CPU). Reuses ecspr_community + the driver's
shared-context solver.
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
from run_community_ecspr import load_axes, solve_pairs_on_graph  # noqa: E402

ELEMENTS = ["C", "N", "S", "P"]


def bh_q(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def _g(reff):
    return 0.0 if (not np.isfinite(reff) or reff <= 0) else 1.0 / reff


def run_null(real_w, real_t, null_w, null_t, bipartite_dir, axes, out_tsv):
    real_orgs = sorted(real_w)
    null_orgs = sorted(null_w)
    ordered = [(a, b) for a in real_orgs for b in real_orgs if a != b]
    print(f"[null] real={real_orgs}  null MAGs (n={len(null_orgs)})={null_orgs}")

    rows = []
    for X in ELEMENTS:
        axX = [a for a in axes if a["element"] == X]
        if not axX:
            continue
        G_X = ec.load_bipartite(Path(bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
        atoms = ec.metabolite_atoms(G_X, X)

        def build(members_w, members_t):
            return ec.build_community_graph(X, members_w, members_t, G_X, atoms=atoms)

        # observed pairwise (real A + real B) conductance, both directions
        obs = {}
        for (a, b) in ordered:
            if (b, a) in obs:  # same undirected graph already solved
                pass
            G = build({a: real_w[a], b: real_w[b]},
                      {a: real_t.get(a, {}), b: real_t.get(b, {})})
            eps = []
            for ax in axX:
                eps.append(((ax["axis_id"], a, b), ("met", ax["source"], a), ("met", ax["sink"], b)))
                eps.append(((ax["axis_id"], b, a), ("met", ax["source"], b), ("met", ax["sink"], a)))
            r = solve_pairs_on_graph(G, X, eps)
            for k, v in r.items():
                obs[k] = _g(v)

        # swap configs: every {real member M, null MAG R} graph, both directions
        # store swap_g[(M, R, axis_id, direction)] where direction in {"R->M","M->R"}
        swap_g = {}
        for M in real_orgs:
            for R in null_orgs:
                G = build({M: real_w[M], R: null_w[R]},
                          {M: real_t.get(M, {}), R: null_t.get(R, {})})
                eps = []
                for ax in axX:
                    # null producer -> real consumer  (R -> M)
                    eps.append(((M, R, ax["axis_id"], "R->M"),
                                ("met", ax["source"], R), ("met", ax["sink"], M)))
                    # real producer -> null consumer  (M -> R)
                    eps.append(((M, R, ax["axis_id"], "M->R"),
                                ("met", ax["source"], M), ("met", ax["sink"], R)))
                r = solve_pairs_on_graph(G, X, eps)
                for k, v in r.items():
                    swap_g[k] = _g(v)

        # assemble per (axis, A, B): producer-swap {g(R->B)} + consumer-swap {g(A->R)}
        for ax in axX:
            aid = ax["axis_id"]
            for (A, B) in ordered:
                g_obs = obs.get((aid, A, B), 0.0)
                prod_null = [swap_g.get((B, R, aid, "R->M"), 0.0) for R in null_orgs]
                cons_null = [swap_g.get((A, R, aid, "M->R"), 0.0) for R in null_orgs]
                pooled = np.array(prod_null + cons_null, dtype=float)
                n = len(pooled)
                p = (1 + int(np.sum(pooled >= g_obs))) / (1 + n)
                rows.append(dict(
                    element=X, axis_id=aid, category=ax["category"],
                    source=ax["source"], sink=ax["sink"], producer=A, consumer=B,
                    g_obs=g_obs, null_mean=float(pooled.mean()),
                    null_max=float(pooled.max()), null_frac_ge=float(np.mean(pooled >= g_obs)),
                    p_emp=p))
        print(f"[null] [{X}] {len(axX)} axes x {len(ordered)} pairs done", flush=True)

    df = pd.DataFrame(rows)
    df["q_bh"] = bh_q(df["p_emp"].to_numpy())
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[null] wrote {out_tsv} ({len(df):,} rows)")
    _summary(df)
    return df


def _summary(df):
    sig = df[(df["q_bh"] < 0.1) & (df["g_obs"] > 0)]
    print(f"\n[null] significant handoffs (q<0.1, g_obs>0): {len(sig):,} / {len(df):,}")
    ery = sig[(sig["producer"] == "ERY") | (sig["consumer"] == "ERY")]
    print(f"[null]   involving Erythrobacter: {len(ery):,}")
    print("\n[null] significant handoffs by ordered pair (producer -> consumer):")
    g = (sig.groupby(["producer", "consumer"])
            .agg(n=("axis_id", "nunique"), med_p=("p_emp", "median"))
            .reset_index().sort_values("n", ascending=False))
    for r in g.itertuples(index=False):
        star = "  <-- Erythrobacter" if "ERY" in (r.producer, r.consumer) else ""
        print(f"    {r.producer} -> {r.consumer}: {r.n} sig axes (median p={r.med_p:.3f}){star}")
    print("\n[null] top Erythrobacter-specific handoffs (lowest p, g_obs > null_max):")
    strong = ery[ery["g_obs"] > ery["null_max"]].sort_values("p_emp").head(15)
    for r in strong.itertuples(index=False):
        print(f"    [{r.element}] {r.producer}->{r.consumer} {r.category:10s} "
              f"{r.source}->{r.sink}  g_obs={r.g_obs:.2f} vs null_max={r.null_max:.2f}  "
              f"p={r.p_emp:.3f} q={r.q_bh:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real-weights", required=True)
    ap.add_argument("--real-transporters", required=True)
    ap.add_argument("--null-weights", required=True)
    ap.add_argument("--null-transporters", required=True)
    ap.add_argument("--bipartite-dir", required=True)
    ap.add_argument("--axes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    run_null(pickle.load(open(a.real_weights, "rb")),
             pickle.load(open(a.real_transporters, "rb")),
             pickle.load(open(a.null_weights, "rb")),
             pickle.load(open(a.null_transporters, "rb")),
             Path(a.bipartite_dir), load_axes(Path(a.axes)), Path(a.out))


if __name__ == "__main__":
    main()
