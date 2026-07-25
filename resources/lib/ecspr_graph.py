"""ECSPr as a measuring instrument: an atom-resolved network, two terminals, one solution.

WHAT THIS IS
------------
A graph plus two metabolite terminals in; a :class:`Solution` out. The solution can be
interrogated at any node or metabolite for **current** and for **voltage**. That is the
whole contract. The engine takes no perturbation argument and knows nothing about how the
network was built: a knockout is a network built from a smaller gene set, and the
comparison is a subtraction the caller does. `methods/pulsechase/run_pulsechase.py` has
stated this contract in its own docstring since before it was implementable.

WHY THE ATOM GRAPH AND NOT THE STAR
-----------------------------------
On the incumbent star topology a metabolite is a single node joined to reaction-node hubs.
Eliminating a reaction node -- which is exactly what a Woodbury update does -- leaves a
CLIQUE over every participant, so two participants that share no atom still get a
conductance between them. Measured on MNXR106432 (carbon): the pyruvate-NADPH and
pyruvate-CoA channels, which carry ZERO carbons, each score 0.700 while the real
1-carbon pyruvate-CO2 channel scores 0.033 -- a zero-atom channel 21x more conductive
than a real one. Here a node IS an atom (metabolite, canonical rank) and an edge IS an
atom transfer, so there is no hub to manufacture a channel and an input terminal can carry
an **atom mask** at all.

THE MEASUREMENT THIS EXISTS FOR
-------------------------------
Inject unit current at a growth substrate; merge every biomass precursor into one
aggregate ground; solve once. Then read, per precursor, how much of that current it
actually draws (:meth:`Solution.delivered`). A loss-of-function edge does not merely lower
a point-to-point conductance -- it REDISTRIBUTES flow, and a starved precursor shows up as
its share collapsing while the total barely moves. That question is invisible to a
two-terminal probe, which is why the two-terminal probe was not enough.

Per-precursor readout must be **current**, not potential. Merged precursors all sit at the
same (ground) potential by construction, so no potential-based readout can distinguish
among them; "shifts current away from a precursor" is inherently a current statement.
Voltage is still exposed, because a non-terminal metabolite's atoms sit at genuinely
different potentials and that is worth being able to see.

HOW A TERMINAL WORKS
--------------------
A terminal is a SET of atom nodes, shorted by exact node contraction -- parallel
conductances kept, intra-terminal edges dropped. Not a large-weight approximation.

One deliberate departure from `ecspr_atom_graph.merge_terminals`: parallel edges are NOT
summed during contraction. Every original edge stays its own row and remembers its
ORIGINAL endpoints. That is what makes per-precursor attribution exact -- with the edges
summed there is no way to ask which precursor a merged edge fed -- and it is electrically
identical, since parallel conductances between the same contracted pair carry current in
the same proportion either way.

THE SOLVE
---------
Reused wholesale from :mod:`ecspr_directed`: signed incidence, softplus-smoothed diode,
Newton with an energy-Armijo line search, CHOLMOD symbolic reuse. Nothing about the
rectified network changes -- only the topology it runs on, and the fact that terminals are
node sets rather than single nodes. When ``gm == gp`` everywhere the network is symmetric
and one linear solve is exact; that path is taken explicitly rather than left for Newton
to discover.

Env: numpy + scipy (+ networkx only for the parity referent). No SCADC paths, no canon.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ecspr_directed import (build_incidence, directed_ceff, _SPDReuse, _softmax0,
                            DIODE_BACKWARD_FLOOR, DIODE_SMOOTH_DELTA, DIRECTED_TOL)

# ---------------------------------------------------------------------------
# Tolerances. REFF_EPS / IEFF_EPS / SELFTEST_TOL moved here from the retired
# `ecspr_solver.py`; they governed the star's two-terminal R_eff and now govern
# this lane's. `ecspr_atom_graph.ATOM_REFF_EPS` remains a separate name for the
# separate (undirected, all-paths) atom measurement -- the two are allowed to
# move independently and must not be conflated.
# ---------------------------------------------------------------------------

# Below this a resistance drop is float-cancellation noise, not real reinforcement.
REFF_EPS = 1e-12
IEFF_EPS = 1e-12

# Tolerance the self-tests hold the solver to against an independent dense rebuild.
# The pulse-chase suite's correctness gate imports this bound rather than restating it.
SELFTEST_TOL = 1e-9

# Terminal supernode keys. Real nodes are (metabolite, rank) tuples whose first element is
# an MNXM string, so these can never collide with one.
SRC_SUPERNODE = ("__terminal__", "source")
SNK_SUPERNODE = ("__terminal__", "sink")


# =====================================================================
# The network
# =====================================================================

@dataclass
class AtomGraph:
    """A directed atom-transfer network.

    ``nodes`` are ``(metabolite, canonical_atom_rank)`` keys; ``edges`` is a list of
    ``(tail_index, head_index)`` pairs into ``nodes``; ``gp`` / ``gm`` are the per-edge
    forward / backward conductances (``gm = ratio * gp``, ratio = ``g_rev/g_fwd`` from the
    direction ensemble; ratio 1.0 is the undirected limit).

    Deliberately the same data shape :class:`ecspr_directed.OrientedNet` carries, so it
    hands straight to ``build_incidence`` and ``directed_ceff``. The difference is that
    there are no reaction nodes: an edge IS an atom transfer, so orientation comes from the
    pair table's substrate/product columns and no reac_prop role parsing is needed.

    ``meta`` carries build provenance (element, reaction count, AAM-gap count, ...). The
    builders in :mod:`ecspr_build` fill it; nothing here reads it.
    """
    nodes: list
    edges: list
    gp: np.ndarray
    gm: np.ndarray
    meta: dict = field(default_factory=dict)
    idx: dict = field(default=None, repr=False)
    _by_met: dict = field(default=None, repr=False)

    def __post_init__(self):
        self.gp = np.asarray(self.gp, float)
        self.gm = np.asarray(self.gm, float)
        if self.idx is None:
            self.idx = {nd: i for i, nd in enumerate(self.nodes)}
        if self._by_met is None:
            by = defaultdict(list)
            for nd in self.nodes:
                by[nd[0]].append(nd)
            self._by_met = {m: sorted(v, key=lambda n: n[1]) for m, v in by.items()}

    @property
    def n(self) -> int:
        return len(self.nodes)

    @property
    def m(self) -> int:
        return len(self.edges)

    def metabolites(self) -> list:
        return sorted(self._by_met)

    def atoms_of(self, met: str) -> list:
        """Every atom node of ``met`` present in the graph, ordered by rank. Empty when the
        metabolite is absent -- an absence a caller must be able to see and report, never a
        raise (a precursor missing from a built graph is a coverage fact, not an error)."""
        return list(self._by_met.get(met, ()))

    def ranks_of(self, met: str) -> list:
        return [nd[1] for nd in self._by_met.get(met, ())]

    # -- convenience -------------------------------------------------------
    @classmethod
    def from_edge_records(cls, records, meta=None):
        """Build from ``[(tail_key, head_key, gp, gm), ...]``. Node order is first-seen."""
        nodes, idx, edges, gp, gm = [], {}, [], [], []

        def _i(k):
            j = idx.get(k)
            if j is None:
                j = idx[k] = len(nodes)
                nodes.append(k)
            return j

        for a, b, p, mn in records:
            edges.append((_i(a), _i(b)))
            gp.append(float(p))
            gm.append(float(mn))
        return cls(nodes, edges, np.asarray(gp, float), np.asarray(gm, float),
                   dict(meta or {}), idx)

    # -- serialisation -----------------------------------------------------
    def save(self, path):
        """One ``.npz`` per graph: node keys, edge endpoints, gp, gm, and the meta as JSON.

        Deliberately not a pickle. A staged reference artifact is read by a container that
        may not hold this module, and a pickle would make the artifact depend on the class
        that wrote it; arrays plus JSON depend on nothing.
        """
        import json as _json
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        e = np.asarray(self.edges, np.int64).reshape(-1, 2) if self.edges else \
            np.zeros((0, 2), np.int64)
        meta = {k: v for k, v in self.meta.items() if not hasattr(v, "shape")
                and k != "edge_reactions"}
        np.savez_compressed(
            path,
            node_met=np.array([n[0] for n in self.nodes], dtype=object),
            node_rank=np.array([n[1] for n in self.nodes], dtype=np.int64),
            tail=e[:, 0], head=e[:, 1], gp=self.gp, gm=self.gm,
            meta=np.array(_json.dumps(meta, default=str)))
        return path

    @classmethod
    def load(cls, path):
        import json as _json
        z = np.load(Path(path), allow_pickle=True)
        nodes = list(zip(z["node_met"].tolist(), (int(r) for r in z["node_rank"])))
        edges = list(zip(z["tail"].tolist(), z["head"].tolist()))
        return cls(nodes, edges, z["gp"], z["gm"], _json.loads(str(z["meta"])))

    def with_ratios(self, ratio: float):
        """A copy with every backward ratio forced to ``ratio``. ``ratio=1.0`` is the
        symmetric limit -- the undirected referent the parity gate measures against."""
        return AtomGraph(list(self.nodes), list(self.edges), self.gp.copy(),
                         ratio * self.gp, dict(self.meta), dict(self.idx))


# =====================================================================
# Terminals
# =====================================================================

@dataclass(frozen=True)
class Terminal:
    """A labelled SET of atom nodes, shorted into one super-node by the solve.

    ``Terminal.metabolite(g, mnxm, mask=...)`` is one metabolite's atoms, optionally
    restricted to a subset of canonical ranks -- the **atom mask**, the thing that only
    exists on an atom-resolved graph. ``Terminal.merge(g, [mnxm, ...])`` is the aggregate
    ground: every listed metabolite's atoms in one terminal.
    """
    label: str
    nodes: frozenset
    metabolites: tuple = ()
    missing: tuple = ()          # metabolites named but absent from the graph

    @classmethod
    def metabolite(cls, graph: AtomGraph, mnxm: str, mask=None, label=None) -> "Terminal":
        atoms = graph.atoms_of(mnxm)
        if mask is not None:
            keep = set(int(r) for r in mask)
            atoms = [a for a in atoms if a[1] in keep]
        return cls(label or mnxm, frozenset(atoms), (mnxm,),
                   () if atoms else (mnxm,))

    @classmethod
    def merge(cls, graph: AtomGraph, mnxms, mask=None, label="ground") -> "Terminal":
        """The aggregate ground. ``mask`` is ``{mnxm: [ranks]}`` when given."""
        mnxms = tuple(mnxms)
        atoms, missing = set(), []
        for m in mnxms:
            a = graph.atoms_of(m)
            if mask and m in mask:
                keep = set(int(r) for r in mask[m])
                a = [x for x in a if x[1] in keep]
            if not a:
                missing.append(m)
            atoms.update(a)
        return cls(label, frozenset(atoms), mnxms, tuple(missing))

    @classmethod
    def of_nodes(cls, label: str, nodes) -> "Terminal":
        nodes = frozenset(nodes)
        return cls(label, nodes, tuple(sorted({n[0] for n in nodes})))

    def __len__(self):
        return len(self.nodes)

    def __bool__(self):
        return bool(self.nodes)


# =====================================================================
# The solution
# =====================================================================

class Solution:
    """The result of one solve: potentials, edge currents, and lookups by node or metabolite.

    Sign conventions, stated once because everything downstream depends on them:

    * unit current is injected at the source terminal and drawn at the sink terminal, so
      :attr:`injected` is 1.0 by construction;
    * :meth:`current` returns net **inflow**, so it is ``+injected`` at the sink terminal,
      ``-injected`` at the source terminal, and 0 at every interior node by KCL;
    * :meth:`voltage` is the grounded potential (the sink super-node is the gauge, so its
      potential is 0), and :attr:`total` is ``1 / (V_source - V_sink)`` -- the two-terminal
      effective conductance between the two terminal SETS.

    Every per-node quantity is computed from the ORIGINAL edge endpoints, not the
    contracted ones, which is what lets a precursor merged into the ground still report the
    current it individually draws.
    """

    __slots__ = ("graph", "source", "sink", "total", "injected", "converged",
                 "_phi_c", "_cidx", "_cmap", "_cur", "_otail", "_ohead", "_nin",
                 "_oedge", "iters", "used_newton", "note")

    def __init__(self, graph, source, sink, *, total, phi_c, cidx, cmap, cur,
                 otail, ohead, converged, oedge=None, iters=0, used_newton=True, note=""):
        self.graph = graph
        self.source = source
        self.sink = sink
        self.total = float(total)
        self.injected = 1.0
        self.converged = bool(converged)
        self.iters = int(iters)
        self.used_newton = bool(used_newton)
        self.note = note
        self._phi_c = phi_c            # potentials, indexed by CONTRACTED index
        self._cidx = cidx              # contracted node key -> contracted index
        self._cmap = cmap              # original node key -> contracted node key
        self._cur = cur                # per-edge current, aligned with _otail/_ohead
        self._otail = otail            # original tail node index per kept edge
        self._ohead = ohead            # original head node index per kept edge
        self._oedge = (np.arange(len(cur), dtype=np.int64) if oedge is None
                       else np.asarray(oedge, np.int64))
        self._nin = graph.n

    def original_edge_index(self) -> np.ndarray:
        """Index into ``graph.edges`` for each solved edge -- the trace back from a current
        to the reactions that built its conductance."""
        return self._oedge

    def edge_currents(self) -> tuple:
        """``(original_edge_index, current)`` for every edge that survived contraction."""
        return self._oedge, self._cur

    # -- voltage -----------------------------------------------------------
    def voltage(self, node) -> float:
        """Potential at an original atom node. ``nan`` when the node is outside the
        component the two terminals share -- floating, not zero, and saying so."""
        ck = self._cmap.get(node)
        if ck is None:
            return float("nan")
        j = self._cidx.get(ck)
        return float("nan") if j is None else float(self._phi_c[j])

    def drop(self, a, b) -> float:
        """Potential difference ``V(a) - V(b)`` between two original nodes."""
        return self.voltage(a) - self.voltage(b)

    def voltage_metabolite(self, mnxm) -> dict:
        """Per-atom potentials plus a current-weighted mean.

        A metabolite's atoms sit at genuinely different potentials unless the metabolite is
        inside a terminal, so a bare scalar would be a lie. The weighted mean uses each
        atom's throughput as its weight -- the potential the carbon actually experiences.
        """
        atoms = self.graph.atoms_of(mnxm)
        per = {a[1]: self.voltage(a) for a in atoms}
        wts = {a[1]: self.throughput(a) for a in atoms}
        finite = [(per[r], wts[r]) for r in per if per[r] == per[r]]
        wsum = sum(w for _, w in finite)
        mean = (sum(v * w for v, w in finite) / wsum) if wsum > 0 else (
            float(np.mean([v for v, _ in finite])) if finite else float("nan"))
        return dict(per_atom=per, weighted_mean=float(mean), n_atoms=len(atoms))

    # -- current -----------------------------------------------------------
    def _mask(self, nodes) -> np.ndarray:
        mk = np.zeros(self._nin, dtype=bool)
        gi = self.graph.idx
        for nd in nodes:
            j = gi.get(nd)
            if j is not None:
                mk[j] = True
        return mk

    def _net_inflow(self, nodes) -> float:
        """Net current entering a SET of original nodes.

        Edges with both endpoints inside the set cancel term-by-term, so this is exactly
        the boundary flux. Edges shorted away by contraction (both endpoints inside ONE
        terminal) carry no current in the contracted model and are absent here -- which is
        why the per-precursor currents still sum to the injected current.
        """
        if self._cur.size == 0:
            return 0.0
        mk = self._mask(nodes)
        return float(np.sum(self._cur * (mk[self._ohead].astype(float)
                                         - mk[self._otail].astype(float))))

    def _boundary_flux(self, nodes) -> tuple:
        """``(inflow, outflow)`` across the boundary of a node set, both non-negative."""
        if self._cur.size == 0:
            return 0.0, 0.0
        mk = self._mask(nodes)
        c = self._cur * (mk[self._ohead].astype(float) - mk[self._otail].astype(float))
        return float(np.sum(np.maximum(c, 0.0))), float(np.sum(np.maximum(-c, 0.0)))

    def _boundary_abs(self, nodes) -> float:
        """Current HANDLED by a node set: ``max(inflow, outflow)``.

        Not ``0.5 * sum |i_e|``. The half-sum is right only where inflow equals outflow --
        an interior node -- and is a factor-of-two undercount at a terminal, where all the
        current arrives and none leaves. ``max`` gives the pass-through node its throughput
        and the sink its full draw, which is what "how much carbon goes through here" means
        at both.
        """
        i, o = self._boundary_flux(nodes)
        return max(i, o)

    def current(self, node) -> float:
        """Net current flowing INTO an original atom node (0 at interior nodes by KCL)."""
        return self._net_inflow((node,))

    def throughput(self, node) -> float:
        """Current handled by a node: ``max(inflow, outflow)`` over its incident edges.

        The quantity to read at a CENTRAL metabolite, where the net current is 0 by KCL and
        therefore says nothing -- what matters there is how much carbon goes through, not
        how much stops.
        """
        return self._boundary_abs((node,))

    # -- metabolite-level wrappers ----------------------------------------
    def delivered(self, mnxm) -> float:
        """Net current drawn by a metabolite -- the load-bearing readout.

        For a precursor merged into the aggregate ground this is exactly "how much of the
        injected current this precursor draws". Summed over the ground's metabolites it
        equals :attr:`injected` to solver tolerance (the conservation self-test).
        """
        return self._net_inflow(self.graph.atoms_of(mnxm))

    def throughput_metabolite(self, mnxm) -> float:
        return self._boundary_abs(self.graph.atoms_of(mnxm))

    def share(self, mnxm) -> float:
        """:meth:`delivered` as a fraction of the injected current."""
        return self.delivered(mnxm) / self.injected if self.injected else 0.0

    def delivered_all(self, mnxms) -> dict:
        return {m: self.delivered(m) for m in mnxms}

    # -- diagnostics -------------------------------------------------------
    def conservation_error(self) -> float:
        """``|sum_over_ground_metabolites(delivered) - injected|``."""
        tot = sum(self.delivered(m) for m in self.sink.metabolites)
        return abs(tot - self.injected)

    def __repr__(self):
        return (f"Solution(total={self.total:.6g}, {self.source.label}->{self.sink.label}, "
                f"converged={self.converged})")


def _zero_solution(graph, source, sink, note):
    empty_i = np.zeros(0, dtype=np.int64)
    return Solution(graph, source, sink, total=0.0, phi_c=np.zeros(0),
                    cidx={}, cmap={}, cur=np.zeros(0), otail=empty_i, ohead=empty_i,
                    converged=True, used_newton=False, note=note)


# =====================================================================
# solve
# =====================================================================

def _contract(graph: AtomGraph, source: Terminal, sink: Terminal):
    """Exact node contraction of the two terminals. Parallel edges are NOT summed.

    Returns ``(cnodes, cidx, cmap, ctail, chead, gp, gm, otail, ohead)`` where the ``c*``
    arrays index the contracted node list and the ``o*`` arrays index the ORIGINAL one.
    """
    cmap = {}
    for nd in source.nodes:
        cmap[nd] = SRC_SUPERNODE
    for nd in sink.nodes:
        cmap[nd] = SNK_SUPERNODE

    cnodes, cidx = [SRC_SUPERNODE, SNK_SUPERNODE], {SRC_SUPERNODE: 0, SNK_SUPERNODE: 1}

    def _ci(key):
        j = cidx.get(key)
        if j is None:
            j = cidx[key] = len(cnodes)
            cnodes.append(key)
        return j

    nodes = graph.nodes
    ct, ch, gp, gm, ot, oh, oe = [], [], [], [], [], [], []
    for e, (a, b) in enumerate(graph.edges):
        ka, kb = nodes[a], nodes[b]
        ca, cb = cmap.get(ka, ka), cmap.get(kb, kb)
        if ca == cb:                       # intra-terminal / self-loop: shorted away
            continue
        ct.append(_ci(ca)); ch.append(_ci(cb))
        gp.append(graph.gp[e]); gm.append(graph.gm[e])
        ot.append(a); oh.append(b); oe.append(e)
    # every original node needs a contracted image for the voltage lookup
    for nd in nodes:
        if nd not in cmap:
            cmap[nd] = nd
    return (cnodes, cidx, cmap, np.asarray(ct, np.int64), np.asarray(ch, np.int64),
            np.asarray(gp, float), np.asarray(gm, float),
            np.asarray(ot, np.int64), np.asarray(oh, np.int64),
            np.asarray(oe, np.int64))


def _component(n, tail, head, seed):
    """Node indices reachable from ``seed`` over the (undirected) contracted edges.

    Reachability is asked from the SOURCE super-node only -- seeding it with both terminals
    would make "the sink is reachable" true by construction and silently turn a
    disconnected pair into a finite conductance.
    """
    A = sp.coo_matrix((np.ones(tail.size), (tail, head)), shape=(n, n)).tocsr()
    lab = sp.csgraph.connected_components(A, directed=False)[1]
    return np.flatnonzero(lab == lab[seed])


def _undirected_phi(B, gp, s, t, g):
    """Grounded potential of the plain undirected network -- exact when ``gm == gp``, and
    the warm start otherwise (a warm start changes only the Newton iteration count, never
    the answer: the smoothed energy is strictly convex with a unique grounded minimiser)."""
    n = B.shape[1]
    keep = np.arange(n) != g
    L = (B.T @ sp.diags(gp) @ B).tocsc()
    rhs = np.zeros(n); rhs[s] = 1.0; rhs[t] = -1.0
    phi = np.zeros(n)
    phi[keep] = splu(L[keep][:, keep].tocsc(), permc_spec="MMD_AT_PLUS_A").solve(rhs[keep])
    return phi


def solve(graph: AtomGraph, source: Terminal, sink: Terminal, *, tol=None, warm=True,
          floor=DIODE_BACKWARD_FLOOR, delta=DIODE_SMOOTH_DELTA) -> Solution:
    """Measure ``graph`` between two terminals. No perturbation argument, by design.

    The terminals are shorted by exact contraction, unit current is injected source->sink,
    and the rectified network is solved by the smoothed-diode Newton in
    :mod:`ecspr_directed`. When every backward conductance equals its forward conductance
    the network is symmetric and one linear solve is exact -- that path is taken explicitly.

    An empty terminal, or two terminals with no connecting path, gives a DEFINITE ZERO
    total rather than a raise: the atom graph is disconnected in general and an
    unreachable precursor is a coverage fact the caller must be able to report.
    """
    if source.nodes & sink.nodes:
        raise ValueError(
            f"terminals overlap on {len(source.nodes & sink.nodes)} atom node(s): "
            f"{source.label} and {sink.label} would short into one node")
    if not source.nodes or not sink.nodes:
        return _zero_solution(graph, source, sink, "empty terminal")

    (cnodes, cidx, cmap, ct, ch, gp, gm, ot, oh, oe) = _contract(graph, source, sink)
    if ct.size == 0:
        return _zero_solution(graph, source, sink, "no edges after contraction")

    sub = _component(len(cnodes), ct, ch, 0)
    inc = np.zeros(len(cnodes), bool)
    inc[sub] = True
    if not inc[1]:
        return _zero_solution(graph, source, sink, "terminals disconnected")

    # Restrict to the shared component: the grounded Hessian is singular otherwise.
    keepe = inc[ct]
    pos = np.full(len(cnodes), -1, np.int64)
    pos[sub] = np.arange(sub.size)
    remap = {int(u): int(pos[u]) for u in sub}
    tail = pos[ct[keepe]]
    head = pos[ch[keepe]]
    gp_k, gm_k = gp[keepe], gm[keepe]
    ot_k, oh_k, oe_k = ot[keepe], oh[keepe], oe[keepe]
    si, ti = remap[0], remap[1]

    B = build_incidence(list(zip(tail.tolist(), head.tolist())), len(sub))
    # Ground at the SINK super-node: a gauge choice, and the one that guarantees the
    # grounded node is in the solved component (directed_ceff otherwise grounds index 0
    # unconditionally, which on a contracted graph need not be connected to the terminals).
    gm_eff = np.maximum(gm_k, floor * gp_k)
    symmetric = bool(np.allclose(gm_eff, gp_k, rtol=1e-12, atol=0.0))

    if symmetric:
        phi = _undirected_phi(B, gp_k, si, ti, ti)
        converged, iters, used_newton = True, 0, False
    else:
        phi0 = _undirected_phi(B, gp_k, si, ti, ti) if warm else None
        kw = dict(g=ti, floor=floor, delta=delta, phi0=phi0, return_phi=True,
                  return_iters=True, return_converged=True,
                  reuse=_SPDReuse(B, np.arange(len(sub)) != ti))
        if tol is not None:
            kw["tol"] = tol
        _, iters, phi, converged = directed_ceff(B, gp_k, gm_k, si, ti, **kw)
        used_newton = True

    dv = float(phi[si] - phi[ti])
    total = 1.0 / dv if abs(dv) > REFF_EPS else 0.0

    # Edge currents from the SMOOTHED law -- the same one the solve minimises. Using the
    # ideal diode instead makes KCL fail against the converged phi at ~1e-4 and reads like
    # a solver bug. At gm == gp the smoothing is an exact no-op (i_e = gp * x).
    x = B @ phi
    hx = _softmax0(x, delta)
    cur = gp_k * hx + gm_eff * (x - hx)

    # Potentials, indexed by contracted key, restricted to the solved component.
    cidx_solved = {cnodes[u]: remap[u] for u in sub}
    return Solution(graph, source, sink, total=total, phi_c=np.asarray(phi, float),
                    cidx=cidx_solved, cmap=cmap, cur=np.asarray(cur, float),
                    otail=ot_k, ohead=oh_k, oedge=oe_k, converged=converged, iters=iters,
                    used_newton=used_newton)


def derive_ieff(r_base: float, r_aug: float) -> tuple:
    """(delta_ieff, g_base, g_aug) exact column transform, no re-solve.

    Moved verbatim from the retired ``ecspr_solver``; the incumbent significance chain and
    the pulse-chase suite both read it.
    """
    rb = max(r_base, IEFF_EPS)
    ra = max(r_aug, IEFF_EPS)
    g_base = 1.0 / rb
    g_aug = 1.0 / ra
    d = g_aug - g_base
    if d < IEFF_EPS:
        d = 0.0
    return d, g_base, g_aug


def _reff_dense(G, s, t, wk=None) -> float:
    """Two-terminal effective resistance by a dense grounded solve on a networkx graph.

    Moved verbatim from the retired ``ecspr_solver``. Shares no line of reasoning with the
    Newton/Woodbury paths, which is exactly what makes it the referent the self-tests and
    the pulse-chase correctness gate check against.
    """
    nodes = list(G.nodes); idx = {n: i for i, n in enumerate(nodes)}; n = len(nodes)
    L = np.zeros((n, n))
    for u, v, d in G.edges(data=True):
        w = 1.0 if wk is None else float(d.get(wk, 0.0))
        if w <= 0:
            continue
        i, j = idx[u], idx[v]
        L[i, i] += w; L[j, j] += w; L[i, j] -= w; L[j, i] -= w
    b = np.zeros(n - 1)
    if idx[s] < n - 1:
        b[idx[s]] += 1.0
    if idx[t] < n - 1:
        b[idx[t]] -= 1.0
    phi = np.zeros(n)
    phi[:-1] = np.linalg.solve(L[:-1, :-1], b)
    return float(phi[idx[s]] - phi[idx[t]])


# =====================================================================
# Self-tests
# =====================================================================

def _toy_graph(ratio=1.0):
    """A small atom network: two metabolites feeding a hub, two precursor sinks.

    m0 (2 atoms) -> m1 (2 atoms) -> {p1, p2}, with an asymmetric split so the two
    precursors draw genuinely different currents, plus a dead-end metabolite d0 whose
    net current must be zero and a disconnected metabolite z0.
    """
    rec = []

    def e(a, b, g):
        rec.append((a, b, g, ratio * g))

    e(("m0", 0), ("m1", 0), 1.0)
    e(("m0", 1), ("m1", 1), 1.0)
    e(("m1", 0), ("m1", 1), 0.3)
    e(("m1", 0), ("p1", 0), 2.0)
    e(("m1", 1), ("p2", 0), 0.5)
    e(("m1", 1), ("p2", 1), 0.25)
    e(("m1", 0), ("d0", 0), 0.7)          # dead end: no onward path
    e(("z0", 0), ("z0", 1), 1.0)          # disconnected component
    return AtomGraph.from_edge_records(rec, dict(kind="toy"))


def _selftest_conservation():
    print("[graph] conservation: per-precursor delivered currents sum to the injected current")
    worst = 0.0
    for ratio in (1.0, 0.01):
        g = _toy_graph(ratio)
        src = Terminal.metabolite(g, "m0")
        gnd = Terminal.merge(g, ["p1", "p2"])
        sol = solve(g, src, gnd)
        d1, d2 = sol.delivered("p1"), sol.delivered("p2")
        err = abs(d1 + d2 - sol.injected)
        worst = max(worst, err)
        src_err = abs(sol._net_inflow(src.nodes) + sol.injected)
        dead = abs(sol.delivered("d0"))
        print(f"  ratio={ratio}: p1={d1:.10f} p2={d2:.10f} sum-1={d1+d2-1:.2e} "
              f"src={src_err:.2e} dead-end={dead:.2e} total={sol.total:.6f}")
        assert err < 1e-8, f"conservation broken at ratio={ratio}: {err:.2e}"
        assert src_err < 1e-8, f"source injection wrong at ratio={ratio}: {src_err:.2e}"
        # The dead end's net current IS the KCL residual at an interior node, so it is
        # bounded by the Newton gradient tolerance rather than by machine epsilon: exactly
        # 0 on the symmetric (single linear solve) path, ~1e-8 on the iterated one.
        assert dead < 1e-6, f"dead-end draws current at ratio={ratio}: {dead:.2e}"
        assert d1 > d2 > 0, "the asymmetric split must actually split asymmetrically"
    print(f"  PASS -- worst |sum - injected| = {worst:.2e}\n")
    return worst


def _selftest_single_sink_reduction():
    """A sink terminal of ONE metabolite must reproduce the two-terminal conductance an
    independent dense solve gives on the contracted graph."""
    import networkx as nx
    print("[graph] single-metabolite sink == independent dense two-terminal solve")
    g = _toy_graph(1.0)
    src = Terminal.metabolite(g, "m0")
    for snk_met in ("p1", "p2"):
        snk = Terminal.metabolite(g, snk_met)
        sol = solve(g, src, snk)
        # dense referent on the SAME contraction, summing parallel conductances
        G = nx.Graph()
        for e, (a, b) in enumerate(g.edges):
            ka, kb = g.nodes[a], g.nodes[b]
            ka = "S*" if ka in src.nodes else ka
            kb = "S*" if kb in src.nodes else kb
            ka = "T*" if ka in snk.nodes else ka
            kb = "T*" if kb in snk.nodes else kb
            if ka == kb:
                continue
            if G.has_edge(ka, kb):
                G[ka][kb]["w"] += float(g.gp[e])
            else:
                G.add_edge(ka, kb, w=float(g.gp[e]))
        comp = nx.node_connected_component(G, "S*")
        r = _reff_dense(G.subgraph(comp), "S*", "T*", wk="w")
        err = abs(sol.total - 1.0 / r)
        print(f"  sink={snk_met}: total={sol.total:.10f} dense={1.0/r:.10f} |d|={err:.2e}")
        assert err < SELFTEST_TOL, f"single-sink reduction mismatch |d|={err:.2e}"
    print(f"  PASS (tol {SELFTEST_TOL:.0e})\n")
    return 0.0


def _selftest_symmetric_limit():
    """Ratio forced to 1.0 must reproduce ``ecspr_atom_graph.supernode_ieff`` -- the bridge
    saying this lane is the old undirected atom measurement plus direction."""
    from ecspr_atom_graph import supernode_ieff
    print("[graph] symmetric limit (ratio=1.0) == ecspr_atom_graph.supernode_ieff")
    g = _toy_graph(1.0)
    edges = defaultdict(float)
    for e, (a, b) in enumerate(g.edges):
        u, v = g.nodes[a], g.nodes[b]
        edges[(u, v) if u < v else (v, u)] += float(g.gp[e])
    worst = 0.0
    cases = [(["m0"], ["p1"]), (["m0"], ["p1", "p2"]), (["m0"], ["p2"])]
    for smets, tmets in cases:
        src = Terminal.merge(g, smets, label="S")
        snk = Terminal.merge(g, tmets, label="T")
        mine = solve(g, src, snk).total
        theirs = supernode_ieff(dict(edges), list(src.nodes), list(snk.nodes))
        err = abs(mine - theirs)
        worst = max(worst, err)
        print(f"  {'+'.join(smets)} -> {'+'.join(tmets)}: new={mine:.10f} "
              f"atom_graph={theirs:.10f} |d|={err:.2e}")
        assert err < SELFTEST_TOL, f"symmetric-limit mismatch |d|={err:.2e}"
    print(f"  PASS -- worst |d| = {worst:.2e}\n")
    return worst


def _selftest_atom_mask():
    """A full mask is a no-op; a strict subset can only lower the total."""
    print("[graph] atom mask: full == unmasked, strict subset <= unmasked")
    g = _toy_graph(1.0)
    gnd = Terminal.merge(g, ["p1", "p2"])
    full = solve(g, Terminal.metabolite(g, "m0"), gnd).total
    same = solve(g, Terminal.metabolite(g, "m0", mask=[0, 1]), gnd).total
    part = solve(g, Terminal.metabolite(g, "m0", mask=[0]), gnd).total
    print(f"  unmasked={full:.10f}  mask[0,1]={same:.10f}  mask[0]={part:.10f}")
    assert abs(full - same) < SELFTEST_TOL, "full mask must be a no-op"
    assert part <= full + SELFTEST_TOL, "a strict atom subset cannot conduct better"
    assert part > 0.0, "the masked terminal should still conduct"
    print("  PASS\n")
    return 0.0


def _selftest_disconnected():
    """A terminal with no path to the other returns a definite zero, not a raise."""
    print("[graph] disconnected terminal -> definite zero (no raise)")
    g = _toy_graph(1.0)
    sol = solve(g, Terminal.metabolite(g, "m0"), Terminal.metabolite(g, "z0"))
    print(f"  m0 -> z0: total={sol.total} note={sol.note!r}")
    assert sol.total == 0.0 and sol.converged
    missing = solve(g, Terminal.metabolite(g, "m0"), Terminal.metabolite(g, "nope"))
    print(f"  m0 -> absent metabolite: total={missing.total} note={missing.note!r}")
    assert missing.total == 0.0
    print("  PASS\n")
    return 0.0


def _selftest_voltage():
    """Voltage lookups: the sink terminal is the gauge (0), the source sits at 1/total,
    and an interior metabolite's atoms sit strictly between."""
    print("[graph] voltage: gauge at the sink, source at 1/total, interior in between")
    g = _toy_graph(1.0)
    src = Terminal.metabolite(g, "m0")
    gnd = Terminal.merge(g, ["p1", "p2"])
    sol = solve(g, src, gnd)
    vs = sol.voltage(("m0", 0))
    vt = sol.voltage(("p1", 0))
    vm = sol.voltage_metabolite("m1")
    print(f"  V(source)={vs:.10f}  1/total={1.0/sol.total:.10f}  V(ground)={vt:.10f}")
    print(f"  m1 per-atom={ {k: round(v, 6) for k, v in vm['per_atom'].items()} } "
          f"weighted_mean={vm['weighted_mean']:.6f}")
    assert abs(vt) < 1e-12, "the sink super-node is the gauge"
    assert abs(vs - 1.0 / sol.total) < SELFTEST_TOL, "V(source) must equal 1/total"
    assert abs(sol.drop(("m0", 0), ("p1", 0)) - vs) < SELFTEST_TOL
    for r, v in vm["per_atom"].items():
        assert vt < v < vs, f"interior atom {r} outside the terminal potentials"
    # every atom of a merged terminal is at the same potential -- that is what shorting means
    assert abs(sol.voltage(("p2", 0)) - sol.voltage(("p1", 0))) < 1e-12
    print("  PASS\n")
    return 0.0


def _selftest_throughput():
    """A pass-through metabolite carries the whole injected current; its NET current is 0."""
    print("[graph] throughput at an interior metabolite == injected; net current == 0")
    g = _toy_graph(1.0)
    sol = solve(g, Terminal.metabolite(g, "m0"), Terminal.merge(g, ["p1", "p2"]))
    net = sol.delivered("m1")
    thru = sol.throughput_metabolite("m1")
    print(f"  m1: net={net:.2e}  throughput={thru:.10f}  injected={sol.injected}")
    assert abs(net) < 1e-8, "an interior metabolite must have zero net current (KCL)"
    assert abs(thru - sol.injected) < 1e-6, "everything must pass through m1"
    # at a terminal, throughput is the full draw -- not the half-sum, which undercounts 2x
    for m in ("p1", "p2"):
        assert abs(sol.throughput_metabolite(m) - sol.delivered(m)) < 1e-8, \
            f"{m}: a pure sink's throughput must equal what it draws"
    assert abs(sol.throughput_metabolite("m0") - sol.injected) < 1e-8, \
        "the source's throughput must equal the injected current"
    print("  PASS\n")
    return 0.0


def _selftest_direction_matters():
    """The instrument must actually see direction: a throttled route delivers less."""
    print("[graph] direction: throttling the p1 branch shifts current to p2")
    g = _toy_graph(1.0)
    sym = solve(g, Terminal.metabolite(g, "m0"), Terminal.merge(g, ["p1", "p2"]))
    lof = AtomGraph(list(g.nodes), list(g.edges), g.gp.copy(), g.gm.copy(), dict(g.meta))
    for e, (a, b) in enumerate(lof.edges):
        if lof.nodes[b][0] == "p1":
            lof.gp[e] *= 1e-4                 # a loss-of-function edge, caller-side
            lof.gm[e] = lof.gp[e]
    ko = solve(lof, Terminal.metabolite(lof, "m0"), Terminal.merge(lof, ["p1", "p2"]))
    print(f"  base:  p1={sym.share('p1'):.6f} p2={sym.share('p2'):.6f} total={sym.total:.6f}")
    print(f"  LOF:   p1={ko.share('p1'):.6f} p2={ko.share('p2'):.6f} total={ko.total:.6f}")
    print(f"  d(share p1)={ko.share('p1') - sym.share('p1'):+.6f}  "
          f"d(share p2)={ko.share('p2') - sym.share('p2'):+.6f}")
    assert ko.share("p1") < sym.share("p1") - 0.1, "the LOF must starve p1"
    assert ko.share("p2") > sym.share("p2") + 0.1, "and the carbon must go to p2"
    assert abs(ko.share("p1") + ko.share("p2") - 1.0) < 1e-8
    print("  PASS\n")
    return 0.0


def _selftest() -> int:
    print("=" * 70)
    print("ecspr_graph self-tests")
    print("=" * 70)
    _selftest_conservation()
    _selftest_single_sink_reduction()
    _selftest_symmetric_limit()
    _selftest_atom_mask()
    _selftest_disconnected()
    _selftest_voltage()
    _selftest_throughput()
    _selftest_direction_matters()
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(_selftest())
