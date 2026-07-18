"""ECSPr directed (rectified "diode" network) solver -- the directional companion to
the undirected effective-conductance solver in ``ecspr_solver.py``.

Each reaction edge carries a forward conductance ``g+`` and a backward conductance
``g-``. A metabolite potential field ``phi`` induces a per-edge signed potential drop
``x_e = (B phi)_e`` where ``B`` is the signed incidence (+1 at the tail, -1 at the head).
Current through an edge is *rectified*:

    i_e(x) = g+ * max(x, 0) - g- * max(-x, 0)  =  g+ * max(x, 0) + g- * min(x, 0)

so the edge conducts freely in the forward direction (``g+``) and is throttled in
reverse (``g- << g+`` for an irreversible reaction). Injecting a unit current s->t and
imposing the KKT stationarity conditions is equivalent to minimising the convex, C^1,
piecewise-quadratic energy

    E(phi) = 0.5 * sum_e [ g+ (x_e^+)^2 + g- (x_e^-)^2 ]  -  I^T phi,   I = e_s - e_t

whose unique (grounded) minimiser gives the two-terminal effective conductance

    C_eff = 1 / (phi[s] - phi[t]).

**Symmetric-limit parity.** At ``g- = g+`` the two rectified branches sum to the plain
quadratic ``0.5 x^T diag(g) x - I^T phi`` (because ``max(x,0)^2 + min(x,0)^2 == x^2``),
so the directed solve reproduces the ordinary undirected ``R_eff`` exactly. A directed
model therefore *degrades to the undirected one* precisely where direction is unknown
(ratio 1.0), and any divergence from the undirected answer is attributable to a real
forward/backward asymmetry -- never to the machinery.

**Solver.** Semismooth Newton with an Armijo line search as the globalisation. The
generalised Hessian is ``B^T diag(d) B`` with ``d = g+`` where ``x>0`` else ``g-`` -- a
weighted graph Laplacian that changes only on edges whose active branch flips. It
converges in ~5-10 iterations; an ADMM variant was measured 10-100x slower and dropped.

This is a **pure primitive**: signed incidence + per-edge conductances + terminals in,
``C_eff`` out. It holds no SCADC path, reads no GEM, and does not import ``canon``.
Edge orientation and the ``g-/g+`` ratio are the caller's responsibility -- for ECSPr
they come from reaction roles (substrate -> rxn -> product) and the directionality
ensemble's per-reaction ``g_rev/g_fwd = exp(dG'/RT)`` ratio, both supplied by the SCADC
driver, never by this module.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

# Newton convergence tolerance on the reduced-gradient inf-norm, and iteration cap.
DIRECTED_TOL = 1e-10
DIRECTED_MAXIT = 60

# The grounded potential is unique only when ``g-`` is strictly positive everywhere:
# a *perfect* diode (g- = 0) leaves the throttled side's potential free and the reduced
# Hessian singular. The primitive clips ``g-`` up to this fraction of ``g+`` as a safety
# net -- documented rather than silently accepting a singular system. Callers that want
# a sharper diode should lower this knowingly, not rely on zero.
DIODE_BACKWARD_FLOOR = 1e-9


def build_incidence(edges, n):
    """Signed incidence ``B`` (m x n, CSR): row ``e`` has +1 at the tail and -1 at the
    head. ``edges`` is a sequence of ``(tail_idx, head_idx)`` integer pairs; ``n`` is the
    node count. Orientation is exactly as given -- the caller decides which end is the
    tail (for ECSPr: substrate-side is the tail, product-side is the head)."""
    m = len(edges)
    if m == 0:
        return sp.csr_matrix((0, n))
    rows = np.repeat(np.arange(m), 2)
    cols = np.empty(2 * m, dtype=np.int64)
    vals = np.empty(2 * m, dtype=float)
    for e, (a, b) in enumerate(edges):
        cols[2 * e], cols[2 * e + 1] = a, b
        vals[2 * e], vals[2 * e + 1] = 1.0, -1.0
    return sp.csr_matrix((vals, (rows, cols)), shape=(m, n))


def directed_ceff(B, gp, gm, s, t, g=0, tol=DIRECTED_TOL, maxit=DIRECTED_MAXIT,
                  floor=DIODE_BACKWARD_FLOOR, return_iters=False, phi0=None,
                  return_phi=False):
    """Two-terminal effective conductance of a rectified network by semismooth Newton.

    Parameters
    ----------
    B : (m, n) signed incidence (CSR), from :func:`build_incidence`.
    gp, gm : (m,) forward / backward conductances per edge (``gm`` is clipped up to
        ``floor * gp`` to keep the grounded system nonsingular).
    s, t : source / sink node indices (unit current injected s -> t).
    g : grounded node index (its potential row/col is dropped; default 0).
    phi0 : optional (n,) warm-start potential. The energy is convex with a unique grounded
        minimiser, so ANY start reaches the same answer -- ``phi0`` only changes the
        iteration count, never the result. Seeding from the undirected grounded solution
        (see :meth:`OrientedNet.ceff`) is the pure-scipy stand-in for symbolic factor reuse.

    Returns
    -------
    C_eff (float), or ``(C_eff, iters)`` when ``return_iters`` is True.
    """
    B = B.tocsr() if sp.issparse(B) else sp.csr_matrix(B)
    gp = np.asarray(gp, float)
    gm = np.maximum(np.asarray(gm, float), floor * gp)
    n = B.shape[1]
    I = np.zeros(n); I[s] = 1.0; I[t] = -1.0
    keep = np.arange(n) != g
    phi = np.zeros(n) if phi0 is None else np.array(phi0, float)
    if phi0 is not None:
        phi[g] = 0.0                              # keep the ground pinned

    def cur(x):
        return gp * np.maximum(x, 0.0) + gm * np.minimum(x, 0.0)

    def energy(p):
        x = B @ p
        return 0.5 * np.sum(gp * np.maximum(x, 0.0) ** 2
                            + gm * np.minimum(x, 0.0) ** 2) - I @ p

    nit = 0
    for nit in range(1, maxit + 1):
        x = B @ phi
        grad = B.T @ cur(x) - I
        if np.abs(grad[keep]).max() < tol:
            break
        d = np.where(x > 0.0, gp, gm)
        H = (B.T @ sp.diags(d) @ B).tocsc()[keep][:, keep]
        dphi = np.zeros(n)
        dphi[keep] = splu(H).solve(-grad[keep])
        E0 = energy(phi)
        slope = grad @ dphi
        step = 1.0
        for _ in range(30):                       # Armijo line search = globalisation
            if energy(phi + step * dphi) <= E0 + 1e-4 * step * slope:
                break
            step *= 0.5
        phi = phi + step * dphi

    ceff = 1.0 / (phi[s] - phi[t])
    if return_phi:
        return (ceff, nit, phi) if return_iters else (ceff, phi)
    return (ceff, nit) if return_iters else ceff


# =====================================================================
# Path-free orientation + element solve (lifted from the SCADC driver's
# DirectedNet._build so the engine, not the experiment, owns the primitive).
# Still canon-free and path-free: the caller supplies the graph OBJECT and
# the roles/ratios DICTS; this module reads no path and imports no canon.
# =====================================================================
from dataclasses import dataclass, field  # noqa: E402


@dataclass
class OrientedNet:
    """An oriented rectified network: signed incidence + per-edge fwd/bwd conductances,
    ready to solve any number of (source, sink) pairs. Built by :func:`orient_and_weight`.

    Warm-start reuse: the undirected grounded Laplacian ``L = B^T diag(gp) B`` is factored
    once (lazily, ``splu``) and cached; every directed solve on this net is seeded from the
    undirected potential. Because the rectified energy is convex with a unique grounded
    minimiser, the seed changes only the Newton iteration count, never the answer -- this is
    the pure-scipy stand-in for symbolic factor reuse (no CHOLMOD dependency).
    """
    nodes: list
    idx: dict
    B: object
    gp: np.ndarray
    gm: np.ndarray
    edges: list = field(default_factory=list, repr=False)
    floor: float = DIODE_BACKWARD_FLOOR
    stats: dict = field(default_factory=dict)
    _fac: object = field(default=None, repr=False)

    @property
    def n(self):
        return self.B.shape[1]

    def _undirected_factor(self):
        if self._fac is None:
            L = (self.B.T @ sp.diags(self.gp) @ self.B).tocsc()
            keep = np.arange(self.n) != 0                 # ground node 0, as directed_ceff
            self._fac = (splu(L[keep][:, keep]), keep)
        return self._fac

    def warm_phi(self, s, t):
        """Undirected grounded potential for terminals (s, t) -- the directed warm start."""
        fac, keep = self._undirected_factor()
        I = np.zeros(self.n); I[s] = 1.0; I[t] = -1.0
        phi = np.zeros(self.n)
        phi[keep] = fac.solve(I[keep])
        return phi

    def _terminal(self, k):
        return k if isinstance(k, (int, np.integer)) else self.idx[k]

    def ceff(self, s, t, *, tol=None, warm=True, return_phi=False, return_iters=False):
        si, ti = self._terminal(s), self._terminal(t)
        phi0 = self.warm_phi(si, ti) if warm else None
        kw = dict(floor=self.floor, phi0=phi0, return_phi=return_phi,
                  return_iters=return_iters)
        if tol is not None:
            kw["tol"] = tol
        return directed_ceff(self.B, self.gp, self.gm, si, ti, **kw)

    def reff_undirected(self, s, t, wk=None):
        """Undirected R_eff on the same graph (the incumbent quantity), for parity."""
        from ecspr_solver import _reff_dense
        return _reff_dense(_net_to_nx(self), s if not isinstance(s, int) else self.nodes[s],
                           t if not isinstance(t, int) else self.nodes[t],
                           wk=wk or "w")


def orient_and_weight(G, X, roles, ratios, *, floor=DIODE_BACKWARD_FLOOR,
                      force_noop=False):
    """Orient a met/rxn bipartite base graph into a rectified network (path-free).

    Substrate edge met->rxn, product edge rxn->met (forward current substrate->rxn->product,
    matching the sign of the ensemble's dG'); ``gp = w`` (element weight ``w_{X}``),
    ``gm = ratio * w``. Role-unknown / no-evidence edges stay symmetric (gm=gp) -- a provable
    no-op reproducing the undirected R_eff exactly there. ``force_noop`` forces every ratio
    to 1.0 (the whole net degrades to undirected).

    ``G`` is a networkx-like bipartite graph (``.nodes``, ``.edges(data=True)``); reaction
    nodes are tuples whose ``[0]`` starts with ``"rxn"``. ``roles`` is MNXR -> (S, P) and
    ``ratios`` is MNXR -> g_rev/g_fwd; both are the caller's, supplied as dicts.
    """
    nodes = list(G.nodes)
    idx = {nd: i for i, nd in enumerate(nodes)}
    wk = f"w_{X}"
    edges, gp, gm = [], [], []
    n_dir = n_role_unknown = n_noev = 0
    _is_rxn = lambda nd: isinstance(nd, tuple) and nd and str(nd[0]).startswith("rxn")
    for u, v, d in G.edges(data=True):
        w = float(d.get(wk, 0.0))
        r, m = (u, v) if _is_rxn(u) else (v, u)
        S, P = roles.get(r[1], (set(), set()))
        if m[1] in S:                             # substrate: met -> rxn
            tail, head = idx[m], idx[r]; oriented = True
        elif m[1] in P:                           # product: rxn -> met
            tail, head = idx[r], idx[m]; oriented = True
        else:                                     # role unknown -> symmetric
            tail, head = idx[u], idx[v]; oriented = False; n_role_unknown += 1
        ratio = 1.0 if (force_noop or not oriented) else ratios.get(r[1], 1.0)
        if oriented and not force_noop:
            if r[1] in ratios and ratio != 1.0:
                n_dir += 1
            elif r[1] not in ratios:
                n_noev += 1
        edges.append((tail, head)); gp.append(w); gm.append(ratio * w)
    return OrientedNet(nodes, idx, build_incidence(edges, len(nodes)),
                       np.asarray(gp, float), np.asarray(gm, float), edges, floor,
                       dict(n_edges=len(edges), n_directed_edges=n_dir,
                            n_role_unknown=n_role_unknown, n_no_evidence_edges=n_noev))


def build_ar2m_directed(weights, G_X, base_nodes, wkey, roles, ratios,
                        *, reinforce=True, force_noop=False):
    """Directed analog of ``ecspr_solver.build_ar2m``: ``{mnxr: E_r}`` -> a flat list of
    directed reinforcement edges ``[(tail_key, head_key, gp, gm), ...]`` for :func:`augment`.

    Node identity matches the undirected ``build_ar2m`` byte-for-byte: a reaction the host
    already carries enters as its parallel copy ``("rxn_reinf", mnxr)`` (so conductances add,
    realising ``E_epi300 + E_fosmid`` electrically) when ``reinforce``; a reaction absent
    from the base enters as the canonical ``("rxn", mnxr)``. Topology is always read from
    the canonical node in ``G_X``; edges are kept only to metabolites already in the base --
    exactly the undirected build's rules.

    Each met edge of conductance ``c = w_X * E_r`` is ORIENTED by the reaction's MetaNetX
    roles, the same convention :func:`orient_and_weight` uses on the base: a substrate feeds
    the reaction (``met -> rxn``), a product is fed by it (``rxn -> met``), and
    ``gm = ratio * gp`` throttles the reverse branch (``ratio = g_rev/g_fwd`` from the
    ensemble; ``1.0`` = no evidence -> symmetric). Role-unknown metabolites stay symmetric.
    ``force_noop`` forces every ratio to 1.0 (the whole addition degrades to undirected) --
    the addition-side half of the symmetric-limit parity gate.

    ``roles`` is ``mnxr -> (S, P)`` and ``ratios`` is ``mnxr -> g_rev/g_fwd``, the caller's
    dicts; ``base_nodes`` is the set of node keys already in the oriented base
    (``set(onet.nodes)``). Path-free and canon-free.
    """
    edges = []
    for mnxr, er in weights.items():
        if er <= 0:
            continue
        canon = ("rxn", mnxr)
        if canon not in G_X:
            continue
        in_base = canon in base_nodes
        if in_base:
            if not reinforce:
                continue
            node_id = ("rxn_reinf", mnxr)
        else:
            node_id = canon
        S, P = roles.get(mnxr, (set(), set()))
        ratio = 1.0 if force_noop else ratios.get(mnxr, 1.0)
        for _, m, d in G_X.edges(canon, data=True):
            if m[0] != "met" or m not in base_nodes:
                continue
            w = d.get(wkey, 0.0)
            if w <= 0:
                continue
            c = float(w) * float(er)
            if m[1] in S:                         # substrate: met -> rxn
                edges.append((m, node_id, c, ratio * c))
            elif m[1] in P:                       # product: rxn -> met
                edges.append((node_id, m, c, ratio * c))
            else:                                 # role unknown -> symmetric
                edges.append((node_id, m, c, c))
    return edges


def _net_to_nx(onet):
    """Rebuild a plain undirected nx.Graph (edge weight 'w' = gp) for the parity referent."""
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(onet.nodes)
    for (a, b), w in zip(onet.edges, onet.gp):
        G.add_edge(onet.nodes[a], onet.nodes[b], w=float(w))
    return G


def augment(onet, add_edges):
    """Extend ``onet`` with directed addition edges (the directed analog of build_ar2m).

    ``add_edges``: list of ``(tail_key, head_key, gp, gm)``; node keys not already in
    ``onet.idx`` become new nodes appended after the base indices (e.g. reinforcement
    ``("rxn_reinf", mnxr)`` nodes for a GOF/draw addition). Returns
    ``(B_aug, gp_aug, gm_aug, idx_aug, n_added_nodes)`` where ``n_added_nodes`` is the count
    of newly introduced nodes. For an ECSPr addition every edge joins a NEW reaction node to
    an EXISTING base metabolite, so ``n_added_nodes`` equals the number of distinct reactions
    added -- the ``n_added_rxns`` the undirected report carries.
    """
    idx = dict(onet.idx)
    nodes = list(onet.nodes)
    n_base = len(nodes)

    def _i(k):
        if k not in idx:
            idx[k] = len(nodes); nodes.append(k)
        return idx[k]

    extra = [(_i(a), _i(b)) for a, b, _, _ in add_edges]
    all_edges = list(onet.edges) + extra
    gp_aug = np.concatenate([onet.gp, np.array([e[2] for e in add_edges], float)]) \
        if add_edges else onet.gp
    gm_aug = np.concatenate([onet.gm, np.array([e[3] for e in add_edges], float)]) \
        if add_edges else onet.gm
    B_aug = build_incidence(all_edges, len(nodes))
    return B_aug, gp_aug, gm_aug, idx, len(nodes) - n_base


def directed_solve_element(onet, axes, additions=None, *, warm=True, tol=None):
    """Yield both-lane directed rows per (axis, unit), preserving the reff/ieff schema.

    ``axes``: list of ``(axis_id, s_key, t_key)`` metabolite terminals in ``onet``.
    ``additions``: ``None`` -> base only (delta 0); else a dict ``unit -> [(tail, head, gp,
    gm), ...]`` of directed reinforcement edges. The base directed solve (with its potential)
    is computed once per axis and reused to warm-start every augmented solve.

    Each row carries both lanes: ``r_base/r_aug/delta_reff`` (reff) and
    ``g_base/g_aug/delta_ieff`` (ieff), the same column meanings the undirected solve emits,
    with ``g = C_eff`` (directed conductance) and ``r = 1/C_eff``.
    """
    from ecspr_solver import derive_ieff

    base = {}
    for axis_id, s, t in axes:
        c, phi = onet.ceff(s, t, tol=tol, warm=warm, return_phi=True)
        base[axis_id] = (float(c), phi)

    units = list(additions.items()) if additions else [(None, [])]
    for unit, add_edges in units:
        if add_edges:
            B_a, gp_a, gm_a, idx_a, n_added = augment(onet, add_edges)
        for axis_id, s, t in axes:
            g_base, phi_base = base[axis_id]
            r_base = 1.0 / g_base
            if add_edges:
                si, ti = idx_a[s], idx_a[t]
                phi0 = np.zeros(B_a.shape[1]); phi0[:onet.n] = phi_base   # pad warm start
                kw = dict(floor=onet.floor, phi0=(phi0 if warm else None))
                if tol is not None:
                    kw["tol"] = tol
                g_aug = float(directed_ceff(B_a, gp_a, gm_a, si, ti, **kw))
                r_aug = 1.0 / g_aug
            else:
                g_aug, r_aug, n_added = g_base, r_base, 0
            d_ieff, gb, ga = derive_ieff(r_base, r_aug)
            yield dict(axis_id=axis_id, unit=unit, s=s, t=t,
                       r_base=r_base, r_aug=r_aug, delta_reff=r_base - r_aug,
                       g_base=gb, g_aug=ga, delta_ieff=d_ieff, n_added_rxns=n_added)


def _directed_ceff_dense(B, gp, gm, s, t, g=0, floor=DIODE_BACKWARD_FLOOR):
    """Independent correctness referent for the toy self-tests.

    Minimises the *same* convex, C^1 energy with a quasi-Newton method (L-BFGS-B) that
    shares no code path with the semismooth-Newton solve above -- different algorithm
    family, so agreement is a genuine cross-check, not a tautology. Dense; toy nets only.
    """
    from scipy.optimize import minimize
    B = B.toarray() if sp.issparse(B) else np.asarray(B, float)
    gp = np.asarray(gp, float)
    gm = np.maximum(np.asarray(gm, float), floor * gp)
    n = B.shape[1]
    I = np.zeros(n); I[s] = 1.0; I[t] = -1.0
    free = np.arange(n) != g

    def fg(p_free):
        p = np.zeros(n); p[free] = p_free
        x = B @ p
        E = 0.5 * np.sum(gp * np.maximum(x, 0.0) ** 2
                         + gm * np.minimum(x, 0.0) ** 2) - I @ p
        grad = B.T @ (gp * np.maximum(x, 0.0) + gm * np.minimum(x, 0.0)) - I
        return E, grad[free]

    res = minimize(fg, np.zeros(int(free.sum())), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
    p = np.zeros(n); p[free] = res.x
    return 1.0 / (p[s] - p[t])


# =====================================================================
# Self-tests: (1) Newton vs an independent dense convex-QP referent;
#             (2) symmetric-limit parity vs the undirected _reff_dense.
# =====================================================================

def _toy_diode_nets():
    """A handful of small (edges, n, gp, gm, s, t) diode networks for the referent gate.
    Orientation and ratios are deliberately varied so the active set is nontrivial."""
    nets = []
    # 1) 4-node diamond: two parallel routes, one throttled backward.
    edges = [(0, 1), (0, 2), (1, 3), (2, 3)]
    gp = np.array([1.0, 1.0, 1.0, 1.0]); gm = np.array([1.0, 0.02, 0.5, 0.02])
    nets.append((edges, 4, gp, gm, 0, 3))
    # 2) 5-node chain with a back-edge, mixed ratios.
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 1)]
    gp = np.array([2.0, 1.0, 1.5, 1.0, 0.7]); gm = np.array([0.02, 1.0, 0.03, 0.5, 0.01])
    nets.append((edges, 5, gp, gm, 0, 4))
    # 3) 6-node grid-ish, several irreversible.
    edges = [(0, 1), (1, 2), (0, 3), (3, 4), (4, 2), (2, 5), (1, 4)]
    gp = np.array([1.0, 1.2, 0.8, 1.0, 1.0, 0.6, 0.9])
    gm = np.array([0.01, 0.01, 1.0, 0.02, 1.0, 0.02, 0.5])
    nets.append((edges, 6, gp, gm, 0, 5))
    return nets


def _selftest_directed_referent(tol=1e-8):
    print("[directed] Newton vs independent dense L-BFGS-B referent on toy diode nets")
    worst = 0.0
    for k, (edges, n, gp, gm, s, t) in enumerate(_toy_diode_nets()):
        B = build_incidence(edges, n)
        c_newt, its = directed_ceff(B, gp, gm, s, t, return_iters=True)
        c_dense = _directed_ceff_dense(B, gp, gm, s, t)
        err = abs(c_newt - c_dense)
        worst = max(worst, err)
        print(f"  net {k}: C_eff newton={c_newt:.10f} dense={c_dense:.10f} "
              f"|d|={err:.2e}  ({its} Newton its)")
        assert err < tol, f"net {k}: directed C_eff mismatch |d|={err:.2e}"
    print(f"  PASS -- worst |d| = {worst:.2e}\n")
    return worst


def _selftest_directed_parity():
    """At g- = g+ the directed solve must reproduce the undirected R_eff to solver
    tolerance. Referent is the undirected dense solve already gated in ecspr_solver."""
    import networkx as nx
    from ecspr_solver import _reff_dense, SELFTEST_TOL

    print("[directed] symmetric-limit parity (g- = g+) vs undirected _reff_dense")
    G = nx.Graph()
    G.add_nodes_from([("met", f"m{i}") for i in range(6)])
    G.add_nodes_from([("rxn", f"r{i}") for i in range(5)])
    edges = [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (2, 3),
             (3, 3), (3, 4), (4, 4), (4, 5), (4, 0)]
    nx_edges = [(("rxn", f"r{r % 5}"), ("met", f"m{m}")) for r, m in edges]
    for u, v in nx_edges:
        G.add_edge(u, v, w=1.0)

    nodes = list(G.nodes)
    idx = {nd: i for i, nd in enumerate(nodes)}
    inc = [(idx[u], idx[v]) for u, v in nx_edges]     # orientation is arbitrary here
    B = build_incidence(inc, len(nodes))
    w = np.ones(len(inc))

    worst = 0.0
    pairs = [(("met", "m0"), ("met", "m5")),
             (("met", "m1"), ("met", "m4")),
             (("met", "m2"), ("met", "m5"))]
    for s_nd, t_nd in pairs:
        r_undir = _reff_dense(G, s_nd, t_nd, wk="w")
        c_undir = 1.0 / r_undir
        c_dir = directed_ceff(B, w, w.copy(), idx[s_nd], idx[t_nd])   # g- = g+
        err = abs(c_dir - c_undir)
        worst = max(worst, err)
        print(f"  {s_nd[1]}->{t_nd[1]}: C_dir={c_dir:.10f} C_undir={c_undir:.10f} "
              f"|d|={err:.2e}")
        assert err < SELFTEST_TOL, f"parity mismatch {s_nd}->{t_nd}: |d|={err:.2e}"
    print(f"  PASS -- worst |d| = {worst:.2e} (tol {SELFTEST_TOL:.0e})\n")
    return worst


def _selftest_directed_behaviour():
    """A directed answer must actually differ by direction: an irreversible chain
    conducts much better forward than backward."""
    print("[directed] behavioural: irreversible chain conducts forward >> reverse")
    # 0 ->1 ->2, both edges strongly forward (g- = 0.01 g+).
    edges = [(0, 1), (1, 2)]
    gp = np.array([1.0, 1.0]); gm = np.array([0.01, 0.01])
    B = build_incidence(edges, 3)
    c_fwd = directed_ceff(B, gp, gm, 0, 2)     # inject with the chain
    c_rev = directed_ceff(B, gp, gm, 2, 0)     # inject against the chain
    print(f"  C_fwd={c_fwd:.6f}  C_rev={c_rev:.6f}  ratio={c_fwd / c_rev:.2f}")
    assert c_fwd > 10.0 * c_rev, "forward should conduct >> reverse"
    # At g- = g+ the same chain is symmetric.
    c_fs = directed_ceff(B, gp, gp.copy(), 0, 2)
    c_rs = directed_ceff(B, gp, gp.copy(), 2, 0)
    assert abs(c_fs - c_rs) < 1e-9, "symmetric chain must be direction-blind"
    print(f"  symmetric limit: C_fwd={c_fs:.6f} == C_rev={c_rs:.6f}  PASS\n")
    return c_fwd / c_rev


def _selftest_warm_start(tol=1e-9):
    """A warm start must reach the SAME C_eff as a cold start (convex energy, unique
    grounded minimiser) -- it may only change the iteration count."""
    print("[directed] warm start == cold (phi0 changes only iteration count)")
    worst = 0.0
    for k, (edges, n, gp, gm, s, t) in enumerate(_toy_diode_nets()):
        B = build_incidence(edges, n)
        c_cold, it_cold = directed_ceff(B, gp, gm, s, t, return_iters=True)
        rng = np.random.default_rng(3 + k)
        phi0 = rng.standard_normal(n); phi0[0] = 0.0          # arbitrary grounded start
        c_warm, it_warm = directed_ceff(B, gp, gm, s, t, return_iters=True, phi0=phi0)
        err = abs(c_cold - c_warm); worst = max(worst, err)
        print(f"  net {k}: cold={c_cold:.10f} ({it_cold}it)  warm={c_warm:.10f} "
              f"({it_warm}it)  |d|={err:.2e}")
        assert err < tol, f"net {k}: warm/cold C_eff mismatch |d|={err:.2e}"
    print(f"  PASS -- worst |d| = {worst:.2e}\n")
    return worst


def _selftest_orient_and_solve():
    """orient_and_weight + OrientedNet.ceff + directed_solve_element: at the symmetric
    limit (roles unknown -> gm=gp) the directed solve reproduces the undirected R_eff, and
    the warm-start path (undirected-seeded) hits the same answer."""
    import networkx as nx
    from ecspr_solver import _reff_dense, SELFTEST_TOL

    print("[directed] orient_and_weight/solve_element symmetric-limit == undirected")
    G = nx.Graph()
    G.add_nodes_from([("met", f"m{i}") for i in range(6)])
    G.add_nodes_from([("rxn", f"r{i}") for i in range(5)])
    rm = [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (2, 3),
          (3, 3), (3, 4), (4, 4), (4, 5), (4, 0)]
    for r, m in rm:
        G.add_edge(("rxn", f"r{r % 5}"), ("met", f"m{m}"), w_C=1.0)

    onet = orient_and_weight(G, "C", roles={}, ratios={})     # no roles -> all symmetric
    assert onet.stats["n_role_unknown"] == onet.stats["n_edges"], "expected all-symmetric"
    worst = 0.0
    for s, t in [(("met", "m0"), ("met", "m5")), (("met", "m1"), ("met", "m4"))]:
        c_dir = onet.ceff(s, t)                               # warm-started (undirected seed)
        c_und = 1.0 / _reff_dense(G, s, t, wk="w_C")
        err = abs(c_dir - c_und); worst = max(worst, err)
        print(f"  {s[1]}->{t[1]}: C_dir={c_dir:.10f} C_und={c_und:.10f} |d|={err:.2e}")
        assert err < SELFTEST_TOL, f"{s}->{t}: symmetric-limit mismatch |d|={err:.2e}"

    # directed_solve_element base lane must equal the undirected R_eff at the symmetric limit
    rows = list(directed_solve_element(onet, [("ax", ("met", "m0"), ("met", "m5"))]))
    r_und = _reff_dense(G, ("met", "m0"), ("met", "m5"), wk="w_C")
    err = abs(rows[0]["r_base"] - r_und); worst = max(worst, err)
    assert err < SELFTEST_TOL, f"solve_element r_base parity |d|={err:.2e}"
    # a symmetric addition (gm=gp) must not break the solve or the schema
    add = {"u1": [(("met", "m0"), ("rxn_reinf", "x"), 1.0, 1.0),
                  (("rxn_reinf", "x"), ("met", "m3"), 1.0, 1.0)]}
    arow = list(directed_solve_element(onet, [("ax", ("met", "m0"), ("met", "m5"))], add))[0]
    assert arow["n_added_rxns"] == 1 and np.isfinite(arow["delta_reff"]), "augment broke"
    print(f"  solve_element base parity |d|={err:.2e}; augment delta_reff="
          f"{arow['delta_reff']:.4e}  PASS\n")
    return worst


def _selftest_build_ar2m_directed():
    """``build_ar2m_directed`` at ``force_noop`` must reproduce the undirected
    ``ecspr_solver.build_ar2m`` node identity, connectivity, and conductances exactly --
    the addition-side half of the symmetric-limit invariant. Then the full oriented +
    augmented directed solve at the symmetric limit must equal the undirected R_eff of the
    same augmented graph (base + reinforcement) computed by the incumbent ``_reff_dense``."""
    import networkx as nx
    from ecspr_solver import build_ar2m, _reff_dense, SELFTEST_TOL

    print("[directed] build_ar2m_directed force_noop == undirected build_ar2m + solve parity")
    # A small universe: 3 reactions, one already in the base (-> reinforced parallel copy),
    # one novel, one carrying a role-unknown metabolite; weights on w_C.
    G_X = nx.Graph()
    rxn_mets = {"r0": ["m0", "m1"], "r1": ["m1", "m2"], "r2": ["m2", "m3"]}
    for r, ms in rxn_mets.items():
        for i, m in enumerate(ms):
            G_X.add_edge(("rxn", r), ("met", m), w_C=1.0 + 0.5 * i)
    # Base = r0 only (so r0 reinforces, r1/r2 are novel), on mets m0,m1,m2,m3.
    base = nx.Graph()
    base.add_edge(("rxn", "r0"), ("met", "m0"), w_C=1.0)
    base.add_edge(("rxn", "r0"), ("met", "m1"), w_C=1.5)
    base.add_edge(("rxn", "rb"), ("met", "m1"), w_C=1.0)   # a base-only reaction for LCC
    base.add_edge(("rxn", "rb"), ("met", "m2"), w_C=1.0)
    base.add_edge(("rxn", "rb"), ("met", "m3"), w_C=1.0)
    base_nodes = set(base.nodes)
    weights = {"r0": 2.0, "r1": 1.0, "r2": 3.0}
    roles = {"r0": ({"m0"}, {"m1"}), "r1": ({"m1"}, {"m2"}), "r2": ({"m2"}, {"m3"})}

    und = build_ar2m(weights, G_X, base_nodes, "w_C", reinforce=True)
    und_edges = {(node, m): c for node, mets in und.items() for m, c in mets}
    dz = build_ar2m_directed(weights, G_X, base_nodes, "w_C", roles, {}, force_noop=True)
    dz_edges = {}
    for a, b, gp, gm in dz:
        assert gp == gm, "force_noop must give a symmetric edge"
        node, m = (a, b) if str(a[0]).startswith("rxn") else (b, a)
        dz_edges[(node, m)] = gp
    assert dz_edges.keys() == und_edges.keys(), \
        f"directed/undirected addition connectivity differs: {dz_edges.keys() ^ und_edges.keys()}"
    worst_w = max(abs(dz_edges[k] - und_edges[k]) for k in und_edges)
    assert worst_w < 1e-12, f"addition conductances differ: worst {worst_w:.2e}"
    print(f"  connectivity: {len(und_edges)} addition edges match; worst |dc|={worst_w:.2e}")

    # Full parity: oriented base + directed additions at the symmetric limit == undirected
    # R_eff of base-plus-reinforcement. Build the undirected augmented graph the incumbent
    # way (parallel conductances add) and compare.
    onet = orient_and_weight(base, "C", roles, {}, force_noop=True)
    add = build_ar2m_directed(weights, G_X, base_nodes, "w_C", roles, {}, force_noop=True)
    rows = list(directed_solve_element(onet, [("ax", ("met", "m0"), ("met", "m3"))],
                                       {"u": add}, warm=True))
    r_aug_dir = next(r["r_aug"] for r in rows if r["unit"] == "u")

    Gu = nx.Graph()
    for u, v, d in base.edges(data=True):
        Gu.add_edge(u, v, w=float(d["w_C"]))
    for a, b, gp, gm in add:                          # parallel add -> sum conductances
        node, m = (a, b) if str(a[0]).startswith("rxn") else (b, a)
        if Gu.has_edge(node, m):
            Gu[node][m]["w"] += gp
        else:
            Gu.add_edge(node, m, w=gp)
    r_aug_und = _reff_dense(Gu, ("met", "m0"), ("met", "m3"), wk="w")
    err = abs(r_aug_dir - r_aug_und)
    assert err < SELFTEST_TOL, f"augmented symmetric-limit parity |d|={err:.2e}"
    print(f"  augmented solve: r_aug directed={r_aug_dir:.10f} undirected={r_aug_und:.10f} "
          f"|d|={err:.2e}  PASS\n")
    return worst_w


def _selftest():
    print("=" * 68)
    print("ecspr_directed self-tests")
    print("=" * 68)
    _selftest_directed_referent()
    _selftest_directed_parity()
    _selftest_directed_behaviour()
    _selftest_warm_start()
    _selftest_orient_and_solve()
    _selftest_build_ar2m_directed()
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(_selftest())
