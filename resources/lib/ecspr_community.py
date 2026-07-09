"""ECSPr community extension -- compartmented, transporter-gated multi-organism graph.

The stock ECSPr universe graphs are compartment-free and transport-blind: the
MetaNetX @MNXD tags are stripped, so transport reactions collapse to A=A and are
dropped, and every metabolite lives in one shared namespace. Concatenating two
organisms' reaction sets on that shared namespace fuses their cytoplasms -- central
metabolites diffuse across a membrane that is not modeled. This module RE-INTRODUCES
the compartmentalization ECSPr discarded, gated by an explicit transporter annotation,
so cross-organism metabolite handoff can only occur through transporter-mediated
exchange with a shared extracellular pool.

Topology (confirmed model):
  * per-organism cytoplasm: every metabolite node is namespaced by organism,
        ("met", MNXM, org);  reactions ("rxn", MNXR, org) wire only that org's copies.
        => NO cytoplasmic metabolite node is ever shared between two organisms.
  * one shared extracellular pool ENV: ("met", MNXM, "ENV"), materialized only for
        metabolites at least one member transports (ENV is not the universe).
  * transporter-gated crossing: for organism o and each substrate M it transports,
        an edge ("met", M, o) <-> ("met", M, "ENV") with conductance
        atoms_X(M) * transporter_belief.  No transporter -> no ENV edge -> M cannot
        leave o's cytoplasm.  This edge IS the membrane.

A metabolite therefore flows A -> B iff A exports it to ENV AND B imports it from ENV
(the biologically precise reading of the shared-ENV model). The permissive variant
("at least one partner transports"; lifted from the cyanoverse two-MAG interface
model) is available via `permissive_env=True` as a sensitivity upper bound.

The graph produced here is a plain networkx.Graph with the same node/edge conventions
the verified SMW solver (ecspr_solver.py) consumes, so effective resistance /
conductance between any two (compartmented) metabolite endpoints is read off with
SMWSolver/reff_base UNCHANGED. At three organisms the graph is small; we solve it
directly (no low-rank pair-update needed).

Interface is swappable: `build_community_graph` takes the per-organism reaction
weights and transporter->substrate maps as data, and a `currency_set` /
`permissive_env` knob, so a coarser or richer transporter mapping swaps in behind the
same API without touching the solver.

Env: numpy + networkx (CPU). Reuses ecspr_solver (numpy/scipy).
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ecspr_solver import SMWSolver, reff_base  # noqa: E402

ELEMENTS = ["C", "N", "S", "P"]
ENV = "ENV"


# =====================================================================
# universe helpers
# =====================================================================

def load_bipartite(path: Path, element: str) -> nx.Graph:
    """Load the universe element bipartite pkl; drop edges with w_X <= 0."""
    with open(path, "rb") as fh:
        G = pickle.load(fh)
    wk = f"w_{element}"
    drop = [(u, v) for u, v, d in G.edges(data=True) if d.get(wk, 0) <= 0]
    G.remove_edges_from(drop)
    return G


def metabolite_atoms(G_X: nx.Graph, element: str) -> dict:
    """atoms_X[MNXM] = number of element-X atoms in the metabolite, read from the
    universe graph. The edge weight w_X of ("rxn",r)-("met",m) is a property of the
    METABOLITE (stoichiometry-aware atom count), constant across m's incident edges;
    we take the max over incident edges to be robust to any zero-padding."""
    wk = f"w_{element}"
    atoms: dict = {}
    for u, v, d in G_X.edges(data=True):
        w = float(d.get(wk, 0.0))
        if w <= 0:
            continue
        m = u if (isinstance(u, tuple) and u[0] == "met") else v
        if isinstance(m, tuple) and m[0] == "met":
            key = m[1]
            if w > atoms.get(key, 0.0):
                atoms[key] = w
    return atoms


# =====================================================================
# community graph build
# =====================================================================

def build_community_graph(element: str,
                          org_weights: dict,
                          transporters: dict,
                          G_X: nx.Graph,
                          currency_set: set | None = None,
                          permissive_env: bool = False,
                          atoms: dict | None = None) -> nx.Graph:
    """Build the compartmented, transporter-gated element-X community graph.

    Parameters
    ----------
    element      : "C" | "N" | "S" | "P"
    org_weights  : {org: {MNXR: E_r}}  per-organism reaction evidence weights.
    transporters : {org: {MNXM: belief}}  per-organism transporter->substrate belief.
    G_X          : the element-X universe bipartite graph (edges pruned to w_X>0).
    currency_set : if given, restrict ENV to these MNXM (transportable currency).
    permissive_env : ignored here (affects scoring / disjoint baseline, not build).
    atoms        : optional precomputed metabolite_atoms(G_X, element).

    Returns a networkx.Graph with nodes ("rxn", MNXR, org), ("met", MNXM, org),
    ("met", MNXM, ENV) and weighted key f"w_{element}".
    """
    wk = f"w_{element}"
    if atoms is None:
        atoms = metabolite_atoms(G_X, element)
    G = nx.Graph()

    # ---- per-organism cytoplasm subgraphs (namespaced) ----
    for org, weights in org_weights.items():
        for mnxr, er in weights.items():
            er = float(er)
            if er <= 0:
                continue
            canon = ("rxn", mnxr)
            if canon not in G_X:
                continue
            rnode = ("rxn", mnxr, org)
            for _, m, d in G_X.edges(canon, data=True):
                if not (isinstance(m, tuple) and m[0] == "met"):
                    continue
                w = float(d.get(wk, 0.0))
                if w <= 0:
                    continue
                mnode = ("met", m[1], org)
                G.add_edge(rnode, mnode, **{wk: w * er})

    # ---- shared ENV + transporter-gated crossing ----
    for org, subs in transporters.items():
        for mnxm, belief in subs.items():
            belief = float(belief)
            if belief <= 0:
                continue
            if currency_set is not None and mnxm not in currency_set:
                continue
            ax = atoms.get(mnxm, 0.0)
            if ax <= 0:            # metabolite carries no element-X atoms -> no X path
                continue
            org_node = ("met", mnxm, org)
            env_node = ("met", mnxm, ENV)
            G.add_edge(org_node, env_node, **{wk: ax * belief})
    return G


def assert_compartment_isolation(G: nx.Graph) -> dict:
    """Verify no cytoplasmic metabolite node is shared across organisms and that every
    cross-organism edge runs through ENV. Returns a small stats dict; raises on
    violation. This is the central correctness guarantee of the method."""
    # metabolite node = ("met", mnxm, comp). Cross-organism leakage would require a
    # single node reachable as two organisms' cytoplasm -- impossible by construction
    # (comp is in the key), so we instead verify every edge is either intra-compartment
    # or an org<->ENV transport edge; never org_a<->org_b directly.
    orgs = set()
    n_transport = 0
    for u, v, d in G.edges(data=True):
        cu = u[2] if (isinstance(u, tuple) and len(u) == 3) else None
        cv = v[2] if (isinstance(v, tuple) and len(v) == 3) else None
        comps = {c for c in (cu, cv) if c is not None}
        non_env = comps - {ENV}
        orgs |= non_env
        if ENV in comps:
            n_transport += 1
            continue
        # both endpoints in a cytoplasm: they must be the SAME organism
        if len(non_env) > 1:
            raise AssertionError(f"direct cross-organism edge (no ENV): {u} -- {v}")
    orgs.discard(ENV)
    # explicit node-overlap check: same (kind, id) met node must not appear under two orgs
    # is structurally impossible (comp in key) but assert the ENV pool is the only bridge.
    return {"organisms": sorted(orgs), "n_transport_edges": n_transport,
            "n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges()}


# =====================================================================
# effective resistance / conductance between compartmented endpoints
# =====================================================================

def community_reff(G: nx.Graph, src: tuple, snk: tuple, element: str) -> float:
    """R_eff between two nodes on the community graph via the verified SMW solver.
    Returns +inf if src and snk are not in the same connected component (no path ->
    infinite resistance)."""
    wk = f"w_{element}"
    if src not in G or snk not in G:
        return float("inf")
    # quick same-component gate (solver takes the global LCC; a cross-component pair
    # would raise or silently sit outside the LCC).
    try:
        cc = nx.node_connected_component(G, src)
    except KeyError:
        return float("inf")
    if snk not in cc:
        return float("inf")
    try:
        solver = SMWSolver(G, [src], [snk], edge_weight_key=wk)
    except ValueError:
        return float("inf")
    if not solver.pairs:
        return float("inf")
    return float(reff_base(solver))


def handoff_score(G_comm: nx.Graph, G_base: nx.Graph, src: tuple, snk: tuple,
                  element: str) -> dict:
    """Handoff = what the partnership enables between src and snk.

    delta_conductance = G_eff(community) - G_eff(baseline), where baseline is the
    graph with the partner's contribution removed (transport edges cut, or the
    partner organism dropped -- caller supplies G_base). Conductance G = 1/R_eff;
    an unreachable pair has G=0.  A positive delta means the community routes
    element flow between the endpoints better than the baseline alone.
    """
    r_comm = community_reff(G_comm, src, snk, element)
    r_base = community_reff(G_base, src, snk, element)
    g_comm = 0.0 if not np.isfinite(r_comm) or r_comm <= 0 else 1.0 / r_comm
    g_base = 0.0 if not np.isfinite(r_base) or r_base <= 0 else 1.0 / r_base
    return {"r_comm": r_comm, "r_base": r_base, "g_comm": g_comm, "g_base": g_base,
            "delta_g": g_comm - g_base, "delta_reff": r_base - r_comm}


def metabolite_handoff_map(G: nx.Graph, organisms: list) -> list:
    """Interpretable cross-feeding map: for each ordered pair (A,B) and each shared
    metabolite M, a handoff A->B exists iff A EXPORTS M (M has a transport edge to ENV
    on A's side AND M participates in >=1 of A's reactions -- i.e. A actually makes/uses
    it, not a dead-end importer) and B IMPORTS M (transport edge on B's side AND M
    participates in >=1 of B's reactions). Returns [{producer, consumer, mnxm}].

    This is the structural companion to the conductance metric: it names WHICH
    metabolites are handed off, independent of any axis endpoint choice.
    """
    def _rxn_degree(node):
        # number of reaction neighbors of a metabolite node (excludes ENV transport)
        if node not in G:
            return 0
        return sum(1 for nb in G.neighbors(node)
                   if isinstance(nb, tuple) and nb[0] == "rxn")

    # metabolites present in ENV (transported by >=1 member)
    env_mets = {n[1] for n in G.nodes
                if isinstance(n, tuple) and n[0] == "met" and n[2] == ENV}
    out = []
    for M in env_mets:
        exporters = set()
        importers = set()
        for o in organisms:
            mnode = ("met", M, o)
            if G.has_edge(mnode, ("met", M, ENV)) and _rxn_degree(mnode) >= 1:
                exporters.add(o)   # symmetric transport: an exporter is also an importer
                importers.add(o)
        for A in exporters:
            for B in importers:
                if A != B:
                    out.append({"producer": A, "consumer": B, "mnxm": M})
    return out


def drop_transport_edges(G: nx.Graph) -> nx.Graph:
    """Return a copy of G with all org<->ENV transport edges removed (the disjoint
    union of cytoplasms). Used as the no-exchange baseline."""
    H = G.copy()
    drop = [(u, v) for u, v in H.edges()
            if (isinstance(u, tuple) and u[2] == ENV) or (isinstance(v, tuple) and v[2] == ENV)]
    H.remove_edges_from(drop)
    return H


# =====================================================================
# axis parsing (axis id = category__source__sink__element)
# =====================================================================

def parse_axis_id(axis_id: str) -> dict | None:
    parts = axis_id.split("__")
    if len(parts) != 4:
        return None
    cat, src, snk, elem = parts
    return {"category": cat, "source": src, "sink": snk, "element": elem}


# =====================================================================
# self-test -- toy graphs, no external annotation needed
# =====================================================================

def _toy_universe() -> nx.Graph:
    """A tiny universe: reactions r1..r4 over metabolites S, X, Y, P (element C).
    r1: S<->X   r2: X<->P   r3: S<->Y   r4: Y<->P    (atom count w_C = 1 everywhere)
    """
    G = nx.Graph()
    edges = [("r1", "S"), ("r1", "X"), ("r2", "X"), ("r2", "P"),
             ("r3", "S"), ("r3", "Y"), ("r4", "Y"), ("r4", "P")]
    for r, m in edges:
        G.add_edge(("rxn", r), ("met", m), w_C=1.0)
    return G


def _selftest() -> int:
    G_X = _toy_universe()
    atoms = metabolite_atoms(G_X, "C")
    assert atoms == {"S": 1.0, "X": 1.0, "Y": 1.0, "P": 1.0}, atoms

    # Organism A carries r1 (S<->X); organism B carries r2 (X<->P). Neither alone
    # links S to P across the membrane. A exports X; B imports X.
    org_weights = {"A": {"r1": 1.0}, "B": {"r2": 1.0}}
    transporters = {"A": {"X": 1.0}, "B": {"X": 1.0}}
    G = build_community_graph("C", org_weights, transporters, G_X)
    stats = assert_compartment_isolation(G)
    assert stats["organisms"] == ["A", "B"], stats
    print(f"[selftest] build ok: {stats}")

    # (1) each organism's copy of a metabolite is a DISTINCT node, and the ONLY
    #     bridge between two organisms' copies of the same metabolite is ENV:
    #     deleting the ENV pool must disconnect ("met","X","A") from ("met","X","B").
    xa, xb = ("met", "X", "A"), ("met", "X", "B")
    assert xa in G and xb in G and xa != xb
    assert G.has_edge(xa, xb) is False, "direct cross-organism metabolite edge!"
    H = G.copy()
    H.remove_nodes_from([n for n in H.nodes if isinstance(n, tuple) and n[2] == ENV])
    assert not nx.has_path(H, xa, xb), "cytoplasms fuse without ENV -- membrane leak!"
    print("[selftest] (1) organism cytoplasms isolated; only ENV bridges them: PASS")

    # (2) a metabolite with no transporter (P) has no ENV edge
    assert ("met", "P", ENV) not in G, "P should not be in ENV (no transporter)"
    assert ("met", "X", ENV) in G, "X should be in ENV (both transport it)"
    print("[selftest] (2) untransported metabolite has no ENV node: PASS")

    # (3) cross-organism S(A) -> P(B) is finite in the community, infinite in the
    #     disjoint union (transport edges removed)
    src, snk = ("met", "S", "A"), ("met", "P", "B")
    r_comm = community_reff(G, src, snk, "C")
    G_disj = drop_transport_edges(G)
    r_disj = community_reff(G_disj, src, snk, "C")
    print(f"[selftest] (3) R_eff  community={r_comm:.4f}  disjoint={r_disj}")
    assert np.isfinite(r_comm) and r_comm > 0, r_comm
    assert not np.isfinite(r_disj), "disjoint union must not connect S(A)->P(B)"
    print("[selftest] (3) exchange is transporter-gated (disjoint -> inf): PASS")

    # (4) removing B's X-transporter breaks the handoff even with A still exporting
    G_noBimport = build_community_graph("C", org_weights,
                                        {"A": {"X": 1.0}, "B": {}}, G_X)
    r_noB = community_reff(G_noBimport, src, snk, "C")
    print(f"[selftest] (4) R_eff with B's X-importer removed = {r_noB}")
    assert not np.isfinite(r_noB), "no B import -> no handoff"
    print("[selftest] (4) both partners must transport (A export AND B import): PASS")

    # (5) SMW reff_base agrees with a direct dense Laplacian solve on the community
    from ecspr_solver import _reff_dense
    r_dense = _reff_dense(
        nx.Graph((u, v, {"w": d["w_C"]}) for u, v, d in G.edges(data=True)),
        src, snk, wk="w")
    err = abs(r_comm - r_dense)
    print(f"[selftest] (5) SMW={r_comm:.8f} dense={r_dense:.8f} |d|={err:.2e}")
    assert err < 1e-9, f"SMW vs dense mismatch {err:.2e}"
    print("[selftest] (5) SMW vs dense agreement: PASS")

    # (6) handoff score is positive for the true community vs disjoint baseline
    h = handoff_score(G, G_disj, src, snk, "C")
    print(f"[selftest] (6) handoff {h}")
    assert h["delta_g"] > 0 and h["g_base"] == 0.0
    print("[selftest] (6) handoff score positive vs no-exchange baseline: PASS")

    # (7) metabolite-level handoff map names X as the A<->B handoff (symmetric)
    hm = metabolite_handoff_map(G, ["A", "B"])
    pairs = {(h["producer"], h["consumer"], h["mnxm"]) for h in hm}
    assert ("A", "B", "X") in pairs and ("B", "A", "X") in pairs, hm
    # P is not transported, S/Y not shared -> only X handed off
    assert {h["mnxm"] for h in hm} == {"X"}, hm
    print(f"[selftest] (7) metabolite handoff map = {sorted(pairs)}: PASS")

    print("\n[selftest] ALL PASS")
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if cmd == "selftest":
        sys.exit(_selftest())
    print(f"unknown subcommand: {cmd}", file=sys.stderr)
    sys.exit(2)
