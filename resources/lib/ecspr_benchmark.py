"""Solve the ECSPr X/Y benchmark: perturb, measure, emit signed effects.

    solve   one (facet, element) shard -> observations for every condition on it
    merge   shard outputs -> the single observations.tsv the scorer consumes

WHAT IT EMITS, AND WHY IT IS A DIFFERENCE
-----------------------------------------
The scalar is the SIGNED CHANGE IN CONDUCTANCE, `g_perturbed - g_base`, not the
pilot's log ratio.

A perturbation that creates a route has `g_base = 0`, where a log ratio is
undefined and any sentinel value is an arbitrary rank injection -- and creating
a route is precisely what the gain-of-function arm does. The difference is
finite everywhere, is exactly zero where the perturbation is silent, and carries
the right sign for both arms. The scorer reads ORDER ONLY, so this costs nothing
and buys total definedness.

Silent conditions are emitted as ZERO, never dropped. Dropping them would
discard exactly the redundancy signal the essentiality diagnostic looks for: a
knockout that changes nothing is a measurement, not a missing value.

WHY THE PARALLELISM IS PROCESSES AND NOT THREADS
------------------------------------------------
The solve factorizes with SuperLU, which is serial, so thread-level parallelism
inside one solve buys nothing. Worse, leaving the maths libraries unpinned lets
each worker spawn its own thread pool and oversubscribe the node -- measured at
more than ten times slower in this project's null builder. So:

  * `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS` are pinned to 1 BEFORE numpy is
    imported, which is why the pinning sits at the very top of this file and not
    in `main()`. Forked workers inherit the pinned environment.
  * Parallelism comes from a fork-based process pool over CONDITIONS.
  * `--workers 1` is the serial reference. Output is an ordered map, so a
    parallel run is byte-identical to it. That equality is the port's own gate.

THE UNIT OF WORK IS A CONDITION, NOT AN EDGE
--------------------------------------------
Each condition needs one graph rebuild and one sparse factorization, then a
cheap solve per panel edge against that shared factorization. Sharding by edge
would redo the factorization for every cell and is the obvious wrong split.
"""
from __future__ import annotations

import os

# BEFORE numpy. See the module docstring -- this is load-bearing, not tidiness.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import math  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import networkx as nx  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ecspr_solver as _s  # noqa: E402

ELEMENTS = ["C", "N", "S", "P"]

# Module-level worker state. A fork pool inherits this by copy-on-write, so the
# base graph and the universe are read ONCE per shard rather than once per
# condition -- they are tens of megabytes and re-reading them per task would
# dominate the runtime.
_W: dict = {}


def _split(cell) -> list[str]:
    import re
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return []
    return [x.strip() for x in re.split(r"[;,]", str(cell)) if x.strip()]


def base_minus_reactions(base: nx.Graph, mnxrs) -> nx.Graph:
    """Copy of `base` with the given reaction nodes removed -- the graph-level
    effect of a gene knockout."""
    G = base.copy()
    G.remove_nodes_from([("rxn", m) for m in mnxrs if ("rxn", m) in G])
    return G


def base_plus_reactions(base: nx.Graph, universe: nx.Graph, el: str,
                        addmap: dict) -> nx.Graph:
    """Copy of `base` with `addmap = {mnxr: E}` added, weighted w_X(rxn,met)*E
    from the universe. New metabolites ARE introduced -- a heterologous
    substrate absent from the host is the point of a GOF addition.

    A reaction the host ALREADY carries enters as a distinct parallel-copy node
    ("rxn_reinf", mnxr) rather than overwriting the base reaction's edges: extra
    ORF support adds parallel conductance. Collapsing it instead makes the add a
    silent no-op for every in-base reaction, which is a real defect this project
    already hit once -- an entire condition scored as "no effect" because the
    add path quietly did nothing.
    """
    G = base.copy()
    wk = f"w_{el}"
    for mnxr, E in addmap.items():
        canon = ("rxn", mnxr)
        if canon not in universe or E <= 0:
            continue
        node = ("rxn_reinf", mnxr) if canon in base else canon
        for _, m, d in universe.edges(canon, data=True):
            if m[0] != "met":
                continue
            w = d.get(wk, 0.0)
            if w > 0:
                G.add_edge(node, m, **{wk: float(w) * float(E)})
    return G


def conductances(G: nx.Graph, el: str, pairs: list[tuple[str, str]]) -> list[float]:
    """G_eff = 1/R_eff for each (src_mnxm, sink_mnxm), on ONE shared factorization.

    An endpoint outside the largest connected component scores 0.0 -- infinite
    resistance. That is the same value a connected-but-non-conducting pair
    returns, which is exactly why the panel builder gates anchor liveness
    separately: here the two are genuinely indistinguishable and must be.
    """
    wk = f"w_{el}"
    ctx = _s.SMWGraphContext(G, device="cpu", dtype="float64", edge_weight_key=wk)
    lcc = ctx.lcc
    out = []
    for src_m, snk_m in pairs:
        src, snk = ("met", src_m), ("met", snk_m)
        if src not in lcc or snk not in lcc or src == snk:
            out.append(0.0)
            continue
        try:
            solver = _s.SMWSolver(G, [src], [snk], edge_weight_key=wk, context=ctx)
            r = _s.reff_base(solver)
        except (ValueError, ZeroDivisionError, FloatingPointError):
            out.append(0.0)
            continue
        out.append((1.0 / r) if (r and r > 0 and math.isfinite(r)) else 0.0)
    return out


def _init_worker(base_path: str, universe_path: str, el: str,
                 pairs: list, edge_ids: list, e_weight: dict) -> None:
    with open(base_path, "rb") as fh:
        base = pickle.load(fh)
    universe = None
    if universe_path:
        with open(universe_path, "rb") as fh:
            universe = pickle.load(fh)
    _W.update(base=base, universe=universe, el=el, pairs=pairs,
              edge_ids=edge_ids, e_weight=e_weight,
              g_base=conductances(base, el, pairs))


def _one_condition(rec: dict) -> tuple[str, list[float]]:
    base, universe, el = _W["base"], _W["universe"], _W["el"]
    add = {m: _W["e_weight"].get(m, 1.0) for m in rec["add"]}
    G = base
    if rec["dele"]:
        G = base_minus_reactions(G, rec["dele"])
    if add:
        if universe is None:
            # A GOF condition without the universe cannot be applied at all.
            # Returning zeros would report "no effect" for a perturbation that
            # was never performed, which is the worst possible failure here.
            raise RuntimeError(
                f"{rec['condition_id']} inserts {len(add)} reaction(s) but no "
                f"universe graph was supplied -- pass --universe. Emitting "
                f"zeros would score an unapplied perturbation as a null result."
            )
        G = base_plus_reactions(G, universe, el, add)
    if G is base:
        # nothing to do: the perturbation touches no reaction in this element's
        # graph. The effect is exactly zero, and it is EMITTED, not dropped.
        return rec["condition_id"], [0.0] * len(_W["pairs"])
    g = conductances(G, el, _W["pairs"])
    return rec["condition_id"], [b - a for a, b in zip(_W["g_base"], g)]


def cmd_solve(a) -> int:
    key = Path(a.key_dir)
    panel = pd.read_csv(key / "panel_edges.tsv", sep="\t", dtype=str)
    cond = pd.read_csv(key / "conditions.tsv", sep="\t", dtype=str).fillna("")

    panel = panel[panel.element == a.element]
    cond = cond[cond.element == a.element]
    if panel.empty or cond.empty:
        print(f"no panel edges or conditions on element {a.element}; nothing to do")
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["facet", "condition_id", "edge_id", "effect"]) \
          .to_csv(a.out, sep="\t", index=False)
        return 0

    pairs = list(zip(panel.anchor_mnxm, panel.product_mnxm))
    edge_ids = list(panel.edge_id)

    # Per-reaction evidence weight, so a GOF add-back restores a reaction to the
    # conductance the host build actually gives it. Net A is uniform E=1.0, so an
    # absent entry correctly defaults to 1.0 rather than to zero.
    e_weight = {}
    if a.reactions:
        rx = pd.read_csv(a.reactions, sep="\t", dtype=str).fillna("")
        rx = rx[rx.build_id == a.facet] if "build_id" in rx.columns else rx
        e_weight = {r.mnxr: float(r.E) for r in rx.itertuples()
                    if getattr(r, "E", "") not in ("", None)}

    if a.limit:
        # A TEST FLAG, and the output says so. A truncated run that reads like a
        # full one is how a partial result gets reported as a complete
        # measurement.
        cond = cond.head(a.limit)
        print(f"  --limit {a.limit}: TRUNCATED RUN, not a scorable result")

    recs = [dict(condition_id=c.condition_id,
                 add=_split(getattr(c, "add_mnxr", "")),
                 dele=_split(getattr(c, "del_mnxr", "")))
            for c in cond.itertuples()]

    # conditions.tsv does not carry the reaction lists -- they live in the
    # observation tables, keyed by obs_id. Join them here rather than widening
    # the answer key, which must stay a statement of biology.
    if a.observations:
        obs_dir = Path(a.observations)
        amap, dmap = {}, {}
        for f, acol, dcol in (("gof_observations.tsv", "add_mnxr", "del_mnxr"),
                              ("lof_observations.tsv", None, "del_mnxr")):
            p = obs_dir / f
            if not p.exists():
                continue
            o = pd.read_csv(p, sep="\t", dtype=str).fillna("")
            for r in o.itertuples():
                if acol and getattr(r, acol, ""):
                    amap[r.obs_id] = _split(getattr(r, acol))
                if getattr(r, dcol, ""):
                    dmap[r.obs_id] = _split(getattr(r, dcol))
        oid = dict(zip(cond.condition_id, cond.obs_id)) if "obs_id" in cond.columns else {}
        for rec in recs:
            o = oid.get(rec["condition_id"], "")
            rec["add"] = rec["add"] or amap.get(o, [])
            rec["dele"] = rec["dele"] or dmap.get(o, [])

    n_touch = sum(1 for r in recs if r["add"] or r["dele"])
    print(f"{a.facet}/{a.element}: {len(recs)} conditions "
          f"({n_touch} perturb a reaction), {len(edge_ids)} panel edges, "
          f"{a.workers} worker(s)")

    base_path = str(Path(a.base_dir) / a.facet / f"base_{a.element}.pkl")
    uni_path = str(Path(a.universe) / f"mnx_bipartite_{a.element}.pkl") if a.universe else ""

    init_args = (base_path, uni_path, a.element, pairs, edge_ids, e_weight)
    rows = []
    if a.workers <= 1:
        # THE SERIAL REFERENCE. A parallel run must be byte-identical to this.
        _init_worker(*init_args)
        results = [_one_condition(r) for r in recs]
    else:
        with ProcessPoolExecutor(max_workers=a.workers, initializer=_init_worker,
                                 initargs=init_args) as ex:
            # ordered map -- chunks keep each worker on contiguous work while the
            # result order stays the input order, so output does not depend on
            # completion order.
            results = list(ex.map(_one_condition, recs,
                                  chunksize=max(1, len(recs) // (a.workers * 4))))

    for cid, effects in results:
        for eid, eff in zip(edge_ids, effects):
            rows.append((a.facet, cid, eid, eff))

    out = pd.DataFrame(rows, columns=["facet", "condition_id", "edge_id", "effect"])
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    nz = int((out.effect != 0).sum())
    print(f"  wrote {a.out}  {len(out):,} rows, {nz:,} non-zero "
          f"({100 * nz / max(1, len(out)):.1f}%)")
    return 0


def cmd_merge(a) -> int:
    parts = []
    for p in sorted(Path(a.shard_dir).rglob("*.tsv")):
        df = pd.read_csv(p, sep="\t")
        if len(df):
            parts.append(df)
    if not parts:
        print("no shard output found -- refusing to write an empty observations.tsv")
        return 1
    out = pd.concat(parts, ignore_index=True)
    out = out.sort_values(["facet", "condition_id", "edge_id"]).reset_index(drop=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"merged {len(parts)} shard(s) -> {a.out}")
    print(f"  {len(out):,} rows | {out.facet.nunique()} facets | "
          f"{out.condition_id.nunique():,} conditions | {out.edge_id.nunique():,} edges")
    print(f"  non-zero: {int((out.effect != 0).sum()):,} "
          f"({100 * (out.effect != 0).mean():.1f}%)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("solve", help="one (facet, element) shard")
    s.add_argument("--key-dir", required=True, help="the frozen Y")
    s.add_argument("--base-dir", required=True)
    s.add_argument("--universe", default="", help="dir of mnx_bipartite_<el>.pkl")
    s.add_argument("--observations", default="", help="dir of gof_/lof_observations.tsv")
    s.add_argument("--reactions", default="", help="X/reactions.tsv, for E weights")
    s.add_argument("--facet", required=True)
    s.add_argument("--element", required=True, choices=ELEMENTS)
    s.add_argument("--workers", type=int,
                   default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
                   help="worker PROCESSES; 1 is the serial reference")
    s.add_argument("--limit", type=int, default=0, help="first N conditions, for tests")
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_solve)

    m = sub.add_parser("merge", help="shards -> observations.tsv")
    m.add_argument("--shard-dir", required=True)
    m.add_argument("--out", required=True)
    m.set_defaults(fn=cmd_merge)

    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
