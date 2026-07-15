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


# =====================================================================
# Graph construction
# =====================================================================

def load_pairs(path: Path, element: str) -> pd.DataFrame:
    df = pd.read_parquet(path)
    return df[df.element == element]


def atom_edges(pairs: pd.DataFrame, weights: dict):
    """{(node_u, node_v): conductance} for the atom-transfer graph.

    node = (metabolite, canonical atom rank). Conductance is E_r per transiting
    atom; parallel transfers between the same atom pair add, which is the same
    electrical convention the star uses for reinforcement.
    """
    out = defaultdict(float)
    for rec in pairs.itertuples(index=False):
        er = weights.get(rec.mnxr, 0.0)
        if er <= 0:
            continue
        si = [int(v) for v in rec.sub_idx.split(",")]
        pi = [int(v) for v in rec.prod_idx.split(",")]
        for a, b in zip(si, pi):
            u, v = (rec.substrate, a), (rec.product, b)
            if u == v:
                continue
            out[(u, v) if u < v else (v, u)] += float(er)
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
    """

    def __init__(self, edges: dict, name: str = ""):
        g = nx.Graph()
        for (u, v), w in edges.items():
            if w > 0:
                g.add_edge(u, v, w=w)
        if not g.number_of_nodes():
            raise ValueError(f"{name}: empty graph")
        lcc = max(nx.connected_components(g), key=len)
        self.g = g
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
    metabolites already in base". Without it, an addition whose endpoints are all
    outside the base LCC forms a component with no path to ground, and the grounded
    Laplacian is exactly singular -- scipy raises rather than lying, which is the one
    mercy here.

    star : an edge is (new_rxn_node, met) -- keep it iff the met is in the base.
           The reaction node is new by construction and attaches through that met.
    atom : an edge is (met_atom, met_atom) with no new hub to hang things from, so
           BOTH endpoints must already be in the base. An addition can reinforce
           routes among host atoms; it cannot invent a metabolite. That is exactly
           the freedom the star lane has, so the lanes stay comparable.
    """
    out = {}
    for (u, v), w in extra.items():
        if atomic:
            if u in grid.idx and v in grid.idx:
                out[(u, v)] = w
        else:
            iu, iv = u in grid.idx, v in grid.idx
            if iu and iv:
                out[(u, v)] = w
            elif iu != iv:
                # one end is the new reaction node; the other must be a base met
                base_end = u if iu else v
                if base_end[0] == "met":
                    out[(u, v)] = w
    return out


def axis_terminals(grid: Grid, met: str, atomic: bool):
    """The node indices an axis endpoint denotes on this grid."""
    if atomic:
        return sorted(grid.idx[n] for n in grid.lcc if n[0] == met)
    node = ("met", met)
    return [grid.idx[node]] if node in grid.idx else []


def score_grid(grid: Grid, atomic: bool, axes: dict, ax_ids: list,
               fos_extra: dict, label: str):
    """delta_ieff per (fosmid, axis): min-route R_eff, one factorization per fosmid."""
    ends = {}
    for ax_id in ax_ids:
        a = axes[ax_id]
        s = axis_terminals(grid, a["source"][0], atomic)
        t = axis_terminals(grid, a["sink"][0], atomic)
        if s and t and set(s) != set(t):
            ends[ax_id] = (s, t)
    if not ends:
        return {}
    terms = sorted({i for s, t in ends.values() for i in s + t})

    def best(Z):
        out = {}
        for ax_id, (S, T) in ends.items():
            r = min(reff_from_z(Z, terms, a, b) for a in S for b in T if a != b)
            out[ax_id] = r
        return out

    t0 = time.time()
    Zb, _ = grid.zcols(terms)
    base = best(Zb)
    rows = {}
    for k, (contig, extra) in enumerate(fos_extra.items()):
        extra = attachable(extra, grid, atomic)
        if not extra:
            continue
        Za, _ = grid.zcols(terms, extra)
        aug = best(Za)
        for ax_id in ends:
            rb, ra = base[ax_id], aug[ax_id]
            gb = 1.0 / max(rb, 1e-12)
            ga = 1.0 / max(ra, 1e-12)
            d = ga - gb
            rows[(contig, ax_id)] = max(d, 0.0)
        if (k + 1) % 50 == 0:
            print(f"    [{label}] {k+1}/{len(fos_extra)}  ({time.time()-t0:.0f}s)",
                  flush=True)
    print(f"    [{label}] {len(rows):,} cells in {time.time()-t0:.0f}s", flush=True)
    return rows
