"""The atom-transfer graph, and the controlled comparison against the star.

WHAT THIS BUILDS
----------------
Two graphs over the SAME reactions, differing only in topology -- which is the only
way to attribute a difference to topology:

  * `atom`  -- nodes are (metabolite, canonical atom rank); edges are atom transfers.
               NO REACTION NODES. A reaction contributes one edge per mapped atom,
               conductance E_r; parallel transfers add.
  * `star`  -- the incumbent topology: ("rxn", r) joined to each ("met", m) with
               conductance w_X(r,m) * E_r. Restricted to the same reaction set.

The star is the thing under suspicion. A reaction node is a hub, and eliminating it --
which is exactly what the Woodbury update does -- leaves a CLIQUE over every
participant, so every pair of participants gets conductance whether or not one atom
passes between them. Measured on MNXR106432 (acetyl-CoA + CO2 + NADPH = pyruvate +
NADP+ + CoA), carbon: pyruvate-NADPH and pyruvate-CoA each get 0.700 while carrying
ZERO carbons, and the real 1-carbon pyruvate-CO2 channel gets 0.033. A zero-atom
channel 21x more conductive than a real one.

WHY BOTH GRAPHS ARE BUILT FROM THE RESTRICTED REACTION SET
----------------------------------------------------------
Only 37.8% of host reactions with evidence survive the atom-pair filter (the rest are
atom-unbalanced, duplicate-metabolite, or unrankable -- see ecspr_atom_pairs). So the
atom graph is not merely a re-topologised star: it also rests on ~half the evidence.
Comparing a full star against a restricted atom graph would confound TOPOLOGY with
EVIDENCE, and the topology question would be unanswerable from the result. Both
graphs therefore use exactly the reactions that have trustworthy pairs.

ENDPOINT SEMANTICS -- the decision the plan flagged as not-to-be-made-silently
------------------------------------------------------------------------------
On the star an axis is one metabolite to another. On the atom graph a metabolite is a
SET of atoms, so "glucose -> acetyl-CoA" has to say WHICH carbons.

Rejected: shorting a terminal's atoms. It is semantically defensible but it makes the
Laplacian axis-dependent, which forfeits the one-factorization-per-fosmid property
that makes any of this affordable.

Rejected: uniform injection over all of a metabolite's atoms. It reintroduces the
very artifact under test -- extracting uniformly from acetyl-CoA's 23 carbons puts
21/23 of the demand on the CoA scaffold, so the number is dominated by the moiety the
carbon never reaches.

Used: R_eff between each source atom and each sink atom, reported as the MINIMUM --
the best available route. That is what "can carbon get from A to B" means, it needs no
arbitrary choice, and it is axis-independent in the Laplacian: every axis is another
RHS column against one factorization. All pairwise resistances among a terminal set T
come from |T| columns, since R_eff(a,b) = Z_aa + Z_bb - 2*Z_ab.

Env: numpy + scipy + pandas + networkx. No SCADC paths.
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
import scipy.sparse as sp
from scipy.sparse.linalg import factorized
import networkx as nx

# Floor below which an atom-lane R_eff is treated as float noise when inverted to a
# conductance. Promoted from a bare `1e-12` literal in `score_grid` so a consumer
# (the pulse-chase suite) can import the bound rather than restate it. This is the
# ATOM-LANE clamp; it is deliberately a SEPARATE name from the star/production
# solver's `ecspr_solver.REFF_EPS`, which governs a different graph. Do not conflate
# the two -- the atom graph and the star graph are different topologies and their
# floors are allowed to move independently.
ATOM_REFF_EPS = 1e-12


# =====================================================================
# Graph construction
# =====================================================================

def load_pairs(path: Path, element: str) -> pd.DataFrame:
    df = pd.read_parquet(path)
    return df[df.element == element]


def _rank_list(v) -> list:
    """A rank field in two schemas: a comma-joined string of ranks (the incumbent extract
    and the toy fixtures pack a whole substrate->product atom list into one row) OR a single
    int (the frozen reference is at atom-pair granularity -- one correspondence per row)."""
    if isinstance(v, str):
        return [int(x) for x in v.split(",")]
    return [int(v)]


def _weight_list(v, n: int) -> list:
    """The per-atom fanout-dilution weight list, matching `_rank_list`'s two schemas: a
    comma-string is split; a scalar (the reference's per-row `pair_w`) is broadcast to n;
    an ABSENT column (None) or a NaN reads as all-1.0 -- so a table with no `pair_w` (the
    incumbent, the star's) is confident-by-default and this stays byte-identical on it."""
    if isinstance(v, str) and v:
        return [float(x) for x in v.split(",")]
    if v is None:
        return [1.0] * n
    try:
        f = float(v)
    except (TypeError, ValueError):
        return [1.0] * n
    if f != f:                       # NaN
        return [1.0] * n
    return [f] * n


def atom_edges(pairs: pd.DataFrame, weights: dict, *, use_confidence: bool = False):
    """{(node_u, node_v): conductance} for the atom-transfer graph.

    node = (metabolite, canonical atom rank). Conductance is E_r * pair_w per transiting
    atom, where `pair_w` is the pair's fanout-dilution weight: 1.0 for a confident mapping
    (an ensemble consensus), 0.5 for a lone-member correspondence, and conf/wsum for a
    disagreement DILUTED across candidates (see `ecspr_atom_pairs` and the AAM combiner).
    That weight IS the ensemble's confidence expressed as conductance -- it is margin-
    conserving (a source atom's candidate weights sum to 1.0, the dilute-not-gap property
    the pulse-chase I1/I3 checks gate), which is exactly why it, not the raw `confidence`
    column, is what folds into the edge. Parallel transfers between the same atom pair add
    (the electrical convention the star uses for reinforcement).

    Two input schemas are accepted transparently (see `_rank_list`/`_weight_list`): the
    incumbent/toy extract (comma-joined ranks, string or absent `pair_w`) and the frozen
    reference (int ranks, scalar `pair_w`). A table without a `pair_w` column reads as
    all-1.0, so this is byte-identical on confident data -- the promotion's no-op.

    `use_confidence` (default OFF) additionally multiplies each row's provenance
    `confidence`. It is a SENSITIVITY knob, not the default weight: for lone-member and
    disagreement rows `confidence` already contains `pair_w`, so multiplying it in
    double-counts the dilution and breaks margin conservation. Left off, conductance is
    the margin-conserving `pair_w`; a table where every `confidence` is 1.0 (or the column
    is absent) makes even the on-path a no-op.
    """
    out = defaultdict(float)
    for rec in pairs.itertuples(index=False):
        er = weights.get(rec.mnxr, 0.0)
        if er <= 0:
            continue
        si = _rank_list(rec.sub_idx)
        pi = _rank_list(rec.prod_idx)
        pw = _weight_list(getattr(rec, "pair_w", None), len(si))
        conf = 1.0
        if use_confidence:
            cf = getattr(rec, "confidence", None)
            if cf is not None:
                f = float(cf)
                if f == f:           # not NaN
                    conf = f
        for a, b, w in zip(si, pi, pw):
            u, v = (rec.substrate, a), (rec.product, b)
            if u == v:
                continue
            out[(u, v) if u < v else (v, u)] += float(er) * float(w) * conf
    return out


def star_edges(G_X: nx.Graph, element: str, rxns, weights: dict):
    """{(node_u, node_v): conductance} for the star, on the SAME reaction set."""
    wk = f"w_{element}"
    out = defaultdict(float)
    for r in rxns:
        rn = ("rxn", r)
        if rn not in G_X:
            continue
        er = weights.get(r, 0.0)
        if er <= 0:
            continue
        for _, m, d in G_X.edges(rn, data=True):
            atom = d.get(wk, 0.0)
            if atom <= 0:
                continue
            out[(rn, m)] += float(atom) * float(er)
    return out


class Grid:
    """A graph as (nodes, index, COO) plus a grounded factorization.

    Deliberately not networkx at solve time: a fosmid's addition only appends edges,
    and the base arrays never change, so they are built once and reused.

    PRUNING. `prune=True` (default) keeps only the largest connected component in the
    factorized arrays (`nodes`/`idx`/`rows`/`cols`/`cond`), because a grounded Laplacian
    over a disconnected graph is singular and the factorization would raise. `prune=False`
    is the no-loss tier: it RETAINS every component, so no reaction's atoms are silently
    dropped -- the property I9 gates. The full graph is always in `self.g` regardless, and
    the all-paths measurement (`supernode_ieff`) reads `self.g`, so the SOLVE never prunes;
    the flag only governs whether the single-factorization `zcols`/`reff_from_z` path is
    available (it needs one connected component). Do not call `zcols` on a `prune=False`
    grid whose graph is disconnected.
    """

    def __init__(self, edges: dict, name: str = "", prune: bool = True):
        g = nx.Graph()
        for (u, v), w in edges.items():
            if w > 0:
                g.add_edge(u, v, w=w)
        if not g.number_of_nodes():
            raise ValueError(f"{name}: empty graph")
        lcc = max(nx.connected_components(g), key=len) if prune else set(g.nodes)
        self.g = g
        self.pruned = prune
        self.lcc = lcc
        self.nodes = sorted(lcc, key=str)
        self.idx = {n: i for i, n in enumerate(self.nodes)}
        self.n = len(self.nodes)
        H = g.subgraph(lcc)
        r, c, w = [], [], []
        for u, v, d in H.edges(data=True):
            r.append(self.idx[u]); c.append(self.idx[v]); w.append(d["w"])
        self.rows = np.asarray(r, dtype=np.int64)
        self.cols = np.asarray(c, dtype=np.int64)
        self.cond = np.asarray(w, dtype=np.float64)
        self.name = name
        self._lu = None

    def _laplacian(self, extra: dict | None = None):
        r, c, w = self.rows, self.cols, self.cond
        n = self.n
        if extra:
            er, ec, ew = [], [], []
            new = {}
            for (u, v), val in extra.items():
                for x in (u, v):
                    if x not in self.idx and x not in new:
                        new[x] = n + len(new)
                iu = self.idx.get(u, new.get(u))
                iv = self.idx.get(v, new.get(v))
                er.append(iu); ec.append(iv); ew.append(val)
            n = n + len(new)
            r = np.concatenate([r, np.asarray(er, dtype=np.int64)])
            c = np.concatenate([c, np.asarray(ec, dtype=np.int64)])
            w = np.concatenate([w, np.asarray(ew, dtype=np.float64)])
        deg = np.zeros(n)
        np.add.at(deg, r, w)
        np.add.at(deg, c, w)
        data = np.concatenate([-w, -w, deg])
        ri = np.concatenate([r, c, np.arange(n)])
        cj = np.concatenate([c, r, np.arange(n)])
        L = sp.coo_matrix((data, (ri, cj)), shape=(n, n)).tocsc()
        return L[1:, 1:].tocsc(), n

    def zcols(self, terms: list, extra: dict | None = None):
        """Z[:, t] for each terminal t. One factorization, |terms| back-solves.

        Returns (Z, m) with Z[i, k] = potential at node i+1 when unit current is
        injected at terms[k] and drawn from the ground. All pairwise R_eff among
        `terms` follows: R(a,b) = Z_aa + Z_bb - 2*Z_ab.
        """
        Lr, m = self._laplacian(extra)
        solve = factorized(Lr)
        B = np.zeros((m - 1, len(terms)))
        for k, t in enumerate(terms):
            if t > 0:
                B[t - 1, k] = 1.0
        return np.asarray(solve(B)), m


def reff_from_z(Z, terms: list, a: int, b: int) -> float:
    """R_eff(a, b) from the terminal columns. Ground (index 0) is implicit."""
    ka, kb = terms.index(a), terms.index(b)
    zaa = Z[a - 1, ka] if a > 0 else 0.0
    zbb = Z[b - 1, kb] if b > 0 else 0.0
    zab = Z[b - 1, ka] if b > 0 else 0.0
    zba = Z[a - 1, kb] if a > 0 else 0.0
    return float(zaa + zbb - zab - zba)


# =====================================================================
# The controlled head-to-head
# =====================================================================

def permute_pairs(pairs: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Shuffle atom maps WITHIN each reaction: same counts, same degree, no pairing.

    The control that decides whether the atom lane measures chemistry or merely AAM
    coverage. It keeps every reaction, every element, and every atom count, and
    destroys only which-atom-goes-where. If effect sizes survive this, the lane is
    reporting that a reaction was mappable -- not what its atoms do. That confound
    already explained ~74% of an earlier element-graph finding, so it is not
    hypothetical.
    """
    rng = np.random.default_rng(seed)
    out = pairs.copy()
    for (mnxr, el), g in pairs.groupby(["mnxr", "element"], sort=False):
        subs = list(g.substrate)
        prods = list(g["product"])
        perm = rng.permutation(len(prods))
        out.loc[g.index, "product"] = [prods[i] for i in perm]
        out.loc[g.index, "prod_idx"] = [list(g.prod_idx)[i] for i in perm]
    return out


def attachable(extra: dict, grid: Grid, atomic: bool) -> dict:
    """Keep only the added edges that actually attach to the host.

    The incumbent does the same thing and says why: build_ar2m keeps "edges only to
    metabolites already in base". Membership is against the WHOLE base graph (`grid.g`),
    not just its LCC: the all-paths measurement (`supernode_ieff`) works on the full
    graph and returns a definite zero for a disconnected endpoint rather than raising, so
    the old 'must be in the LCC or the grounded Laplacian is singular' guard no longer
    applies -- an addition may legitimately reinforce a minority component.

    star : an edge is (new_rxn_node, met) -- keep it iff the met is in the base.
           The reaction node is new by construction and attaches through that met.
    atom : an edge is (met_atom, met_atom) with no new hub to hang things from, so
           BOTH endpoints must already be in the base. An addition can reinforce
           routes among host atoms; it cannot invent a metabolite. That is exactly
           the freedom the star lane has, so the lanes stay comparable.
    """
    base = set(grid.g.nodes)
    out = {}
    for (u, v), w in extra.items():
        if atomic:
            if u in base and v in base:
                out[(u, v)] = w
        else:
            iu, iv = u in base, v in base
            if iu and iv:
                out[(u, v)] = w
            elif iu != iv:
                # one end is the new reaction node; the other must be a base met
                base_end = u if iu else v
                if base_end[0] == "met":
                    out[(u, v)] = w
    return out


def axis_terminals(grid: Grid, met: str, atomic: bool):
    """The node indices an axis endpoint denotes on this grid (LCC-local)."""
    if atomic:
        return sorted(grid.idx[n] for n in grid.lcc if n[0] == met)
    node = ("met", met)
    return [grid.idx[node]] if node in grid.idx else []


# =====================================================================
# The all-paths, edge-preserving measurement -- replaces MIN-over-pairs
# =====================================================================

def _dense_reff(G: nx.Graph, s, t) -> float:
    """Two-terminal effective resistance by a dense grounded solve on G (edge attr 'w').

    G must be connected (the caller restricts to the shared component). Grounds the last
    node and inverts the reduced Laplacian. Shares no line of reasoning with the Woodbury
    /Z path, which is what makes it the referent the pulse-chase gate checks against.
    """
    nodes = list(G.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    L = np.zeros((n, n))
    for u, v, d in G.edges(data=True):
        w = float(d.get("w", 0.0))
        if w <= 0:
            continue
        i, j = idx[u], idx[v]
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    b = np.zeros(n - 1)
    if idx[s] < n - 1:
        b[idx[s]] += 1.0
    if idx[t] < n - 1:
        b[idx[t]] -= 1.0
    phi = np.zeros(n)
    phi[:-1] = np.linalg.solve(L[:-1, :-1], b)
    return float(phi[idx[s]] - phi[idx[t]])


def merge_terminals(edges: dict, groups: list) -> dict:
    """Short each named atom group into one super-node: sum parallel conductances, drop
    intra-group edges. The exact realization of 'short all source atoms, short all sink
    atoms' -- no big-weight approximation."""
    remap = {}
    for name, atoms in groups:
        for a in atoms:
            remap[a] = name
    out = defaultdict(float)
    for (u, v), w in edges.items():
        uu, vv = remap.get(u, u), remap.get(v, v)
        if uu == vv:
            continue
        key = (uu, vv) if str(uu) < str(vv) else (vv, uu)
        out[key] += w
    return dict(out)


def metabolite_atoms(edges: dict, met: str) -> list:
    """Every atom node (met, rank) of `met` present in an edge dict -- ALL of them, across
    components. An axis endpoint is the whole metabolite, so its terminal is its atom SET."""
    seen = set()
    for (u, v) in edges:
        for nnode in (u, v):
            if isinstance(nnode, tuple) and nnode and nnode[0] == met:
                seen.add(nnode)
    return sorted(seen, key=str)


def supernode_ieff(edges: dict, src_atoms, snk_atoms) -> float:
    """The all-paths, edge-preserving effective conductance (Ieff) between two metabolite
    endpoints. Each endpoint is the SET of its atoms, shorted to a super-node, and the
    conductance is the full two-terminal solve -- every path, every edge, combined in
    parallel.

    THIS REPLACES THE MIN-OVER-ATOM-PAIRS SCORER, and the reason is behavioural, not
    cosmetic. MIN reports only the single best (source-atom, sink-atom) route, so it
    neither COMBINES parallel routes -- a second, independent route through different
    atoms leaves it unchanged -- nor is EDGE-PRESERVING: removing a non-best edge that
    still carries source->sink current leaves it unchanged. A pulse-chase measures where
    the whole label can go, not its single fastest atom, so the endpoint must be the whole
    atom set and every route must count. Shorting the endpoint atoms also BRIDGES what
    atom-level components would split apart -- two routes that share no atom node still
    share the source and sink metabolites -- which is what makes a metabolite-to-metabolite
    conductance well defined at all.

    Returns 0.0 (a definite zero-capacity, not an error) when the endpoints have no
    connecting path or when either atom set is empty. Uses the atom-lane conductance floor
    ATOM_REFF_EPS.
    """
    if not src_atoms or not snk_atoms:
        return 0.0
    merged = merge_terminals(edges, [("S*", set(src_atoms)), ("T*", set(snk_atoms))])
    G = nx.Graph()
    for (u, v), w in merged.items():
        if w > 0:
            G.add_edge(u, v, w=w)
    if "S*" not in G or "T*" not in G or not nx.has_path(G, "S*", "T*"):
        return 0.0
    # _dense_reff grounds the last node and inverts the reduced Laplacian, singular on a
    # disconnected graph. Restrict to the component the endpoints share (they do, above).
    comp = nx.node_connected_component(G, "S*")
    H = G.subgraph(comp)
    return 1.0 / max(_dense_reff(H, "S*", "T*"), ATOM_REFF_EPS)


def _endpoint_atoms(node_set: set, met: str, atomic: bool) -> list:
    """The nodes an axis endpoint denotes, across the WHOLE graph. atom: every (met, rank)
    atom node of the metabolite. star: the single ("met", met) node."""
    if atomic:
        return [n for n in node_set
                if isinstance(n, tuple) and len(n) == 2 and n[0] == met
                and isinstance(n[1], (int, np.integer))]
    node = ("met", met)
    return [node] if node in node_set else []


def score_grid(grid: Grid, atomic: bool, axes: dict, ax_ids: list,
               fos_extra: dict, label: str):
    """delta_ieff per (fosmid, axis) via the ALL-PATHS super-node measurement.

    Replaces the MIN-over-atom-pairs scorer with `supernode_ieff` (see its docstring for
    why MIN is behaviourally wrong): each axis endpoint is the whole metabolite's atom
    SET, shorted, and Ieff is the full two-terminal solve on the ENTIRE graph -- every
    route, edge-preserving, and spanning the components atom-level pruning would separate.

    Operates on `grid.g` (the full graph, not the LCC): shorting a metabolite's atoms
    bridges components, and a route that reaches the endpoint through a minority component
    is real capacity the pruned LCC would silently drop. This forfeits the one-
    factorization-per-fosmid property the MIN path had -- a dense solve per (fosmid, axis)
    -- which is the cost of the correctness the pulse-chase II2 gate demands; speeding it
    back up (e.g. a super-node reff read off the shared factorization) is a separate,
    benchmarked step.
    """
    base_edges = {(u, v): d["w"] for u, v, d in grid.g.edges(data=True)}
    base_nodes = set(grid.g.nodes)
    ends = {}
    for ax_id in ax_ids:
        a = axes[ax_id]
        s = _endpoint_atoms(base_nodes, a["source"][0], atomic)
        t = _endpoint_atoms(base_nodes, a["sink"][0], atomic)
        if s and t and set(s) != set(t):
            ends[ax_id] = (a["source"][0], a["sink"][0])
    if not ends:
        return {}

    t0 = time.time()
    base = {ax_id: supernode_ieff(base_edges,
                                  _endpoint_atoms(base_nodes, sm, atomic),
                                  _endpoint_atoms(base_nodes, km, atomic))
            for ax_id, (sm, km) in ends.items()}
    rows = {}
    for k, (contig, extra) in enumerate(fos_extra.items()):
        extra = attachable(extra, grid, atomic)
        if not extra:
            continue
        aug_edges = dict(base_edges)
        for (u, v), w in extra.items():
            key = (u, v) if str(u) < str(v) else (v, u)
            aug_edges[key] = aug_edges.get(key, 0.0) + w
        aug_nodes = base_nodes | {n for e in extra for n in e}
        for ax_id, (sm, km) in ends.items():
            gb = base[ax_id]
            ga = supernode_ieff(aug_edges,
                                _endpoint_atoms(aug_nodes, sm, atomic),
                                _endpoint_atoms(aug_nodes, km, atomic))
            rows[(contig, ax_id)] = max(ga - gb, 0.0)
        if (k + 1) % 50 == 0:
            print(f"    [{label}] {k+1}/{len(fos_extra)}  ({time.time()-t0:.0f}s)",
                  flush=True)
    print(f"    [{label}] {len(rows):,} cells in {time.time()-t0:.0f}s", flush=True)
    return rows
