"""ECSPr current-flow solver -- Effective Conductance with Stoichiometry-aware PRoblem.

Self-contained port of the scadc reaction-network solver into the fabfos library.
It fuses three scadc source modules into one importable + CLI module:

  * `_reference_try1/betweenness/07b_smw_current_flow.py` -- SMWGraphContext /
    SMWSolver (Laplacian build + dense inverse; the Sherman-Morrison-Woodbury base).
  * `04_reaction_network/_smw.py` -- build_ar2m (evidence-weighted addition maps).
  * `04_reaction_network/_reff.py` -- reff_summary / reff_batch (dR_eff = z^T M dz)
    and the ieff column transform (dI_eff = 1/r_aug - 1/r_base).

The betweenness (Delta-b) machinery of 07b is intentionally dropped: ECSPr scores
effective resistance/conductance, not current-flow betweenness. Only the pieces the
reff/ieff solve needs are retained, so the module is CPU-first (GPU optional, same
math). This keeps `_smw.py`'s importlib-by-stem hack out of the picture.

Math (single source->sink pair; reduced Laplacian L, ground row/col 0 dropped;
b = e_s - e_t):
    R_eff_base = b^T L^-1 b = phi0[s] - phi0[t],   phi0 = L^-1 b = Zst[:, pair]
    A fosmid adds reactions -> low-rank update L_aug = L + E Delta E^T. By Woodbury,
    delta_reff = R_eff_base - R_eff_aug = z^T M Delta z >= 0,
      z = E^T phi0 = phi0[touched],  M = (I + Delta E^T L^-1 E)^-1 = (I + Delta Z_ee)^-1.
    delta_ieff = 1/r_aug - 1/r_base (exact column transform; r_aug = r_base - delta_reff).

Env: numpy + scipy + networkx (CPU). torch-CUDA only if device="cuda".
"""
from __future__ import annotations

import sys
from typing import Iterable

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import factorized
import networkx as nx

# Below this a resistance drop is float-cancellation noise, not real reinforcement.
REFF_EPS = 1e-12
IEFF_EPS = 1e-12


# =====================================================================
# SMW base: per-element shared Laplacian + dense inverse (port of 07b)
# =====================================================================

class SMWGraphContext:
    """Per-element shared SMW state (topology + element weights + dense inverse
    Z = L_r^-1). Computed once per element, shared across all axes on the same
    (graph, edge_weight_key). The O(n^3) inverse is the bottleneck, so sharing it
    across an element's axes removes redundant inversions."""

    def __init__(self, G_base: nx.Graph, device: str = "cpu",
                 dtype: str = "float64", edge_weight_key: str | None = None):
        self.device = device
        self.edge_weight_key = edge_weight_key

        if edge_weight_key is not None:
            G_eff = nx.Graph()
            G_eff.add_nodes_from(G_base.nodes(data=True))
            for u, v, d in G_base.edges(data=True):
                w = float(d.get(edge_weight_key, 0))
                if w > 0:
                    G_eff.add_edge(u, v, **{edge_weight_key: w})
        else:
            G_eff = G_base

        lcc = max(nx.connected_components(G_eff), key=len)
        H = G_eff.subgraph(lcc).copy()
        self.lcc = lcc
        # ground = highest-degree node (least likely touched by additions; the
        # two-terminal R_eff is ground-invariant so any choice is correct).
        nodes = sorted(H.nodes, key=lambda n: (-H.degree(n), str(n)))
        self.nodes = nodes
        self.idx = {n: i for i, n in enumerate(nodes)}
        n = len(nodes)
        self.n = n

        rows, cols, conds = [], [], []
        for u, v, d in H.edges(data=True):
            rows.append(self.idx[u])
            cols.append(self.idx[v])
            conds.append(1.0 if edge_weight_key is None
                         else float(d.get(edge_weight_key, 0)))
        self.rows = np.asarray(rows, dtype=np.int64)
        self.cols = np.asarray(cols, dtype=np.int64)
        self.cond = np.asarray(conds, dtype=np.float64)

        deg = np.zeros(n)
        np.add.at(deg, self.rows, self.cond)
        np.add.at(deg, self.cols, self.cond)
        data = np.concatenate([-self.cond, -self.cond, deg])
        ri = np.concatenate([self.rows, self.cols, np.arange(n)])
        cj = np.concatenate([self.cols, self.rows, np.arange(n)])
        L = sp.coo_matrix((data, (ri, cj)), shape=(n, n)).tocsc()
        self.L_r = L[1:, 1:]
        self._solve = factorized(self.L_r)

        self._gpu = (device in ("cuda", "gpu"))
        if self._gpu:
            import torch
            self.torch = torch
            self.tdev = torch.device("cuda")
            self.tdtype = torch.float64 if dtype == "float64" else torch.float32
            L64 = torch.from_numpy(self.L_r.toarray()).to(self.tdev).to(torch.float64)
            Z64 = torch.linalg.inv(L64)
            self.Z = Z64.to(self.tdtype).contiguous()
            del L64, Z64
            torch.cuda.empty_cache()
        else:
            self.torch = None
            self.tdev = None
            self.tdtype = None


class SMWSolver:
    """Per-axis solver: base potential Zst for each (source, sink) pair. Shares an
    SMWGraphContext across an element's axes when one is supplied."""

    def __init__(self, G_base: nx.Graph, sources: Iterable, sinks: Iterable,
                 device: str = "cpu", dtype: str = "float64",
                 edge_weight_key: str | None = None,
                 context: "SMWGraphContext | None" = None):
        if context is not None:
            assert context.edge_weight_key == edge_weight_key, \
                "context's edge_weight_key must match"
            self.device = context.device
            self.edge_weight_key = edge_weight_key
            self.nodes = context.nodes
            self.idx = context.idx
            self.n = context.n
            self.rows = context.rows
            self.cols = context.cols
            self.cond = context.cond
            self.L_r = context.L_r
            self._solve = context._solve
            self._gpu = context._gpu
            if self._gpu:
                self.torch = context.torch
                self.tdev = context.tdev
                self.tdtype = context.tdtype
                self.Z = context.Z
            lcc = context.lcc
            src_in = [s for s in sources if s in lcc]
            snk_in = [t for t in sinks if t in lcc]
            if not src_in or not snk_in:
                raise ValueError("no sources or sinks in LCC of G_base")
            self._init_per_axis(src_in, snk_in)
            return

        ctx = SMWGraphContext(G_base, device=device, dtype=dtype,
                              edge_weight_key=edge_weight_key)
        self.__init__(G_base, sources, sinks, device=device, dtype=dtype,
                      edge_weight_key=edge_weight_key, context=ctx)

    def _init_per_axis(self, src_in: list, snk_in: list) -> None:
        n = self.n
        self.sources_idx = [self.idx[s] for s in src_in]
        self.sinks_idx = [self.idx[t] for t in snk_in]
        self.src_nodes = src_in
        self.snk_nodes = snk_in
        pairs, rhs_cols = [], []
        for si in self.sources_idx:
            for ti in self.sinks_idx:
                if si == ti:
                    continue
                pairs.append((si, ti))
                b = np.zeros(n - 1)
                if si > 0:
                    b[si - 1] += 1.0
                if ti > 0:
                    b[ti - 1] -= 1.0
                rhs_cols.append(b)
        self.pairs = pairs
        if rhs_cols:
            rhs = np.column_stack(rhs_cols)
            try:
                Zst = self._solve(rhs)
            except Exception:
                Zst = np.column_stack([self._solve(c) for c in rhs_cols])
        else:
            Zst = np.zeros((n - 1, 0))
        self.Zst = Zst.astype(np.float64)


# =====================================================================
# Evidence-weighted addition maps (port of _smw.build_ar2m)
# =====================================================================

def build_ar2m(weights: dict, G_X, base_nodes: set, wkey: str,
               reinforce: bool = True) -> dict:
    """{mnxr: E_r} -> {node_id: [(('met', m), w_X*E_r), ...]} the solve consumes.

    A reaction the host already carries enters as a distinct parallel-copy node
    ("rxn_reinf", mnxr) when reinforce=True (parallel conductances add, realizing
    E_epi300 + E_fosmid electrically); topology is always read from the canonical
    ("rxn", mnxr) node in G_X. Edges kept only to metabolites already in base."""
    ar2m: dict = {}
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
        entries = []
        for _, m, d in G_X.edges(canon, data=True):
            if m[0] != "met" or m not in base_nodes:
                continue
            w = d.get(wkey, 0.0)
            if w > 0:
                entries.append((m, float(w) * float(er)))
        if entries:
            ar2m[node_id] = entries
    return ar2m


# =====================================================================
# Effective resistance / conductance (port of _reff)
# =====================================================================

def _norm_ar2m(solver, ar2m: dict):
    """(ar2m_v, touched, tpos): ar2m_v maps rxn -> [(full_idx, c)]; touched is the
    sorted full met indices > 0 (ground dropped); tpos maps full_idx -> position."""
    ar2m_v: dict = {}
    for r, mets in ar2m.items():
        mi_c = []
        for entry in mets:
            if isinstance(entry, tuple) and len(entry) == 2 and isinstance(entry[0], tuple) \
                    and isinstance(entry[1], (int, float)):
                m, c = entry
                c = float(c)
                if c <= 0:
                    continue
            else:
                m, c = entry, 1.0
            if m in solver.idx:
                mi_c.append((solver.idx[m], c))
        if mi_c:
            ar2m_v[r] = mi_c
    tset = set()
    for r, mi_c in ar2m_v.items():
        for i, _ in mi_c:
            if i > 0:
                tset.add(i)
    touched = sorted(tset)
    tpos = {fi: j for j, fi in enumerate(touched)}
    return ar2m_v, touched, tpos


def _build_delta(ar2m_v: dict, tpos: dict, p: int) -> np.ndarray:
    """Delta on the touched space: Delta[i,i]+=c_i; Delta[i,j]-=c_i c_j/d_r."""
    Delta = np.zeros((p, p), dtype=np.float64)
    for r, mi_c in ar2m_v.items():
        d_r = sum(c for _, c in mi_c)
        if d_r <= 0:
            continue
        for mi, c_i in mi_c:
            if mi == 0:
                continue
            Delta[tpos[mi], tpos[mi]] += c_i
        non_ground = [(mi, c) for (mi, c) in mi_c if mi > 0]
        for mi, c_i in non_ground:
            for mj, c_j in non_ground:
                Delta[tpos[mi], tpos[mj]] -= (c_i * c_j) / d_r
    return Delta


def reff_base(solver, pair: int = 0) -> float:
    """R_eff of the base graph for the axis's (s,t) pair = phi0[s] - phi0[t]."""
    if not solver.pairs:
        return 0.0
    si, ti = solver.pairs[pair]
    phi = solver.Zst[:, pair]
    vs = float(phi[si - 1]) if si > 0 else 0.0
    vt = float(phi[ti - 1]) if ti > 0 else 0.0
    return vs - vt


# =====================================================================
# Direct multi-axis solve -- one factorization per fosmid, no Woodbury
# =====================================================================

class DirectSolver:
    """Every axis for one host, by factoring the augmented Laplacian once.

    WHY THIS EXISTS, WHEN THE WOODBURY PATH ALREADY WORKS
    -----------------------------------------------------
    Woodbury is worth its complexity only when `Z = L^-1` is a precomputed DENSE
    inverse and `Z_ee` is therefore a gather. That is true on the GPU path. It is
    NOT true on CPU: `reff_summary` builds `rhs = zeros((n-1, p))` and calls
    `solver._solve(rhs)`, which is **p sparse solves per (axis, fosmid)** with p =
    the number of touched metabolites. The base factorization is shared, but that
    work is redone for every axis even though `Delta` and `Z_ee` do not depend on
    the axis at all -- only `z = Zst[:, pair]` does.

    So on CPU the update math costs more than the problem it is optimizing.
    Measured on carbon (n_LCC=7685, 40 axes, 8 fosmids): Woodbury 88.8s, this 1.36s
    -- 65x, agreeing to 2.8e-13. Extrapolated over the basis: ~45 min -> ~1 min.

    The augmented graph is the base plus one new node per added reaction, each
    joined to the metabolites it touches. That is a BORDERED system: the base block
    is identical for every fosmid, so it is assembled once (from the context's COO
    arrays) and only the border changes. Factor it, then every axis is a back-solve
    against the same LU.

    This is not a numerical approximation of the Woodbury path -- it is the same
    two-terminal R_eff computed without eliminating the reaction nodes first. It
    keeps the reaction node explicit rather than collapsing it into the clique that
    `_build_delta` forms, which is why it also scales to a graph where a dense
    inverse is not merely slow but impossible (atom-resolved carbon, n ~ 6e4, would
    need ~32 GB for Z alone).

    The Woodbury path is deliberately kept: it is independently verified against a
    dense rebuild (~1e-14), so it is a referent for this one, and the two share no
    math.
    """

    def __init__(self, ctx: SMWGraphContext):
        self.ctx = ctx
        self.n = ctx.n
        self.idx = ctx.idx
        self._base = (ctx.rows, ctx.cols, ctx.cond)
        self._lu_base = None

    def _assemble(self, ar2m_v: dict):
        """Augmented Laplacian, grounded at node 0. Base block + border."""
        br, bc, bw = self._base
        n = self.n
        rn = {r: n + i for i, r in enumerate(ar2m_v)}
        m = n + len(rn)
        er, ec, ew = [], [], []
        for r, mi_c in ar2m_v.items():
            ri = rn[r]
            for mi, c in mi_c:
                er.append(mi); ec.append(ri); ew.append(float(c))
        if er:
            R = np.concatenate([br, np.asarray(er, dtype=np.int64)])
            C = np.concatenate([bc, np.asarray(ec, dtype=np.int64)])
            W = np.concatenate([bw, np.asarray(ew, dtype=np.float64)])
        else:
            R, C, W = br, bc, bw
        deg = np.zeros(m)
        np.add.at(deg, R, W)
        np.add.at(deg, C, W)
        data = np.concatenate([-W, -W, deg])
        ri = np.concatenate([R, C, np.arange(m)])
        cj = np.concatenate([C, R, np.arange(m)])
        L = sp.coo_matrix((data, (ri, cj)), shape=(m, m)).tocsc()
        return L[1:, 1:].tocsc(), m

    @staticmethod
    def _rhs(pairs, m):
        B = np.zeros((m - 1, len(pairs)))
        for k, (si, ti) in enumerate(pairs):
            if si > 0:
                B[si - 1, k] += 1.0
            if ti > 0:
                B[ti - 1, k] -= 1.0
        return B

    @staticmethod
    def _harvest(phi, pairs):
        out = []
        for k, (si, ti) in enumerate(pairs):
            vs = float(phi[si - 1, k]) if si > 0 else 0.0
            vt = float(phi[ti - 1, k]) if ti > 0 else 0.0
            out.append(vs - vt)
        return out

    def base_reff(self, pairs: list) -> list:
        """R_eff for every axis on the unaugmented host. Factored once, cached."""
        if self._lu_base is None:
            Lr, _ = self._assemble({})
            self._lu_base = factorized(Lr)
        B = self._rhs(pairs, self.n)
        return self._harvest(np.asarray(self._lu_base(B)), pairs)

    def aug_reff(self, ar2m_v: dict, pairs: list) -> list:
        """R_eff for EVERY axis after one fosmid's additions. One factorization.

        `pairs` are (source, sink) indices into the base node order; the border
        never renumbers the base, so the same indices address the augmented system.
        """
        if not ar2m_v:
            return self.base_reff(pairs)
        Lr, m = self._assemble(ar2m_v)
        solve = factorized(Lr)
        return self._harvest(np.asarray(solve(self._rhs(pairs, m))), pairs)


def reff_direct(dsolver: "DirectSolver", ar2m_v: dict, pairs: list,
                rb: list) -> list:
    """(delta_reff, r_base, r_aug) per axis for one fosmid.

    delta >= 0 by Rayleigh monotonicity -- adding conductance cannot raise a
    two-terminal resistance. It is clamped at REFF_EPS for the same reason the
    Woodbury path is: a delta below float cancellation noise is noise, not
    reinforcement, and letting a -1e-16 through as a negative resistance drop
    would be a lie told to one decimal place.
    """
    ra = dsolver.aug_reff(ar2m_v, pairs)
    out = []
    for r0, r1 in zip(rb, ra):
        d = r0 - r1
        if d < REFF_EPS:
            d = 0.0
        out.append((d, r0, r0 - d))
    return out


def reff_summary(solver, ar2m: dict, pair: int = 0,
                 r_base: float | None = None) -> tuple:
    """(delta_reff, r_base, r_aug, n_added) for one fosmid on one axis (float64)."""
    rb = reff_base(solver, pair) if r_base is None else float(r_base)
    ar2m_v, touched, tpos = _norm_ar2m(solver, ar2m)
    n_added = len(ar2m_v)
    if not touched:
        return 0.0, rb, rb, n_added
    p = len(touched)
    E_red = np.asarray([m - 1 for m in touched], dtype=np.int64)
    Delta = _build_delta(ar2m_v, tpos, p)

    if getattr(solver, "_gpu", False):
        torch = solver.torch
        E_t = torch.as_tensor(E_red, device=solver.tdev)
        Zee = solver.Z.index_select(0, E_t).index_select(1, E_t)
        Z_ee = Zee.detach().to("cpu").numpy().astype(np.float64)
    else:
        rhs = np.zeros((solver.n - 1, p))
        rhs[E_red, np.arange(p)] = 1.0
        try:
            Y_red = solver._solve(rhs)
        except Exception:
            Y_red = np.column_stack([solver._solve(rhs[:, j]) for j in range(p)])
        Z_ee = np.asarray(Y_red)[E_red, :].astype(np.float64)

    z = solver.Zst[E_red, pair].astype(np.float64)
    try:
        M = np.linalg.inv(np.eye(p) + Delta @ Z_ee)
    except np.linalg.LinAlgError:
        M = np.linalg.pinv(np.eye(p) + Delta @ Z_ee)
    delta = float(z @ (M @ (Delta @ z)))
    if delta < REFF_EPS:
        delta = 0.0
    return delta, rb, rb - delta, n_added


def derive_ieff(r_base: float, r_aug: float) -> tuple:
    """(delta_ieff, g_base, g_aug) exact column transform, no re-solve."""
    rb = max(r_base, IEFF_EPS)
    ra = max(r_aug, IEFF_EPS)
    g_base = 1.0 / rb
    g_aug = 1.0 / ra
    d = g_aug - g_base
    if d < IEFF_EPS:
        d = 0.0
    return d, g_base, g_aug


# =====================================================================
# Self-test (port of _reff_selftest.py) -- Woodbury vs dense
# =====================================================================

def _reff_dense(G, s, t, wk=None) -> float:
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


def _selftest() -> int:
    G = nx.Graph()
    G.add_nodes_from([("met", f"m{i}") for i in range(6)])
    G.add_nodes_from([("rxn", f"r{i}") for i in range(5)])
    edges = [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (2, 3),
             (3, 3), (3, 4), (4, 4), (4, 5), (4, 0)]
    for r, m in edges:
        G.add_edge(("rxn", f"r{r % 5}"), ("met", f"m{m}"))
    s, t = ("met", "m0"), ("met", "m5")
    solver = SMWSolver(G, [s], [t], edge_weight_key=None, device="cpu")

    rb_solver = reff_base(solver)
    rb_dense = _reff_dense(G, s, t)
    print(f"R_eff_base: solver={rb_solver:.10f} dense={rb_dense:.10f} "
          f"|d|={abs(rb_solver - rb_dense):.2e}")
    assert abs(rb_solver - rb_dense) < 1e-9, "base R_eff mismatch"

    cases = [
        {("rxn", "new1"): [(("met", "m1"), 1.5), (("met", "m4"), 2.0)]},
        {("rxn", "new2"): [(("met", "m0"), 1.0), (("met", "m5"), 1.0)]},
        {("rxn", "n3"): [(("met", "m2"), 0.7), (("met", "m3"), 0.4), (("met", "m5"), 1.1)]},
        {("rxn", "a"): [(("met", "m1"), 2.0), (("met", "m2"), 2.0)],
         ("rxn", "b"): [(("met", "m3"), 0.5), (("met", "m4"), 0.5)]},
    ]
    worst = 0.0
    for k, ar2m in enumerate(cases):
        d_solver, rb, raug, na = reff_summary(solver, ar2m)
        Ga = G.copy()
        for r, mets in ar2m.items():
            for m, c in mets:
                Ga.add_edge(r, m, w=float(c))
        for u, v, dd in Ga.edges(data=True):
            if "w" not in dd:
                dd["w"] = 1.0
        d_dense = _reff_dense(G, s, t) - _reff_dense(Ga, s, t, wk="w")
        err = abs(d_solver - d_dense)
        worst = max(worst, err)
        di, gb, ga = derive_ieff(rb, raug)
        print(f"case {k}: dR solver={d_solver:.10f} dense={d_dense:.10f} "
              f"|d|={err:.2e} (dI={di:.6f}, n_added={na})")
        assert d_solver >= -1e-12, "delta_reff < 0"
    assert worst < 1e-8, f"delta_reff mismatch (worst {worst:.2e})"
    print(f"\nPASS -- worst |d| = {worst:.2e}")
    return 0


def _selftest_direct() -> int:
    """DirectSolver vs the Woodbury path vs a dense rebuild -- three ways, one answer.

    The three do not share math: Woodbury eliminates the reaction node into a clique
    and updates an inverse; DirectSolver keeps the node and factors the bordered
    system; `_reff_dense` builds the whole thing and solves it. Agreement across all
    three is evidence. Agreement across two implementations of one idea is not.

    Includes a CYCLE among the added edges on purpose: two reactions sharing two
    metabolites close a loop, which is the case that makes an incidence matrix
    rank-deficient and would break a naive edge-space formulation.
    """
    G = nx.Graph()
    G.add_nodes_from([("met", f"m{i}") for i in range(6)])
    G.add_nodes_from([("rxn", f"r{i}") for i in range(5)])
    for r, m in [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (2, 3),
                 (3, 3), (3, 4), (4, 4), (4, 5), (4, 0)]:
        G.add_edge(("rxn", f"r{r % 5}"), ("met", f"m{m}"))
    s, t = ("met", "m0"), ("met", "m5")

    ctx = SMWGraphContext(G, edge_weight_key=None, device="cpu")
    dsolver = DirectSolver(ctx)
    solver = SMWSolver(G, [s], [t], edge_weight_key=None, context=ctx)
    pairs = solver.pairs

    rb_direct = dsolver.base_reff(pairs)
    rb_wood = [reff_base(solver, k) for k in range(len(pairs))]
    rb_dense = _reff_dense(G, s, t)
    print(f"r_base: direct={rb_direct[0]:.12f} woodbury={rb_wood[0]:.12f} "
          f"dense={rb_dense:.12f}")
    assert abs(rb_direct[0] - rb_dense) < 1e-9, "direct base != dense"
    assert abs(rb_wood[0] - rb_dense) < 1e-9, "woodbury base != dense"

    cases = [
        {("rxn", "new1"): [(("met", "m1"), 1.5), (("met", "m4"), 2.0)]},
        {("rxn", "new2"): [(("met", "m0"), 1.0), (("met", "m5"), 1.0)]},
        {("rxn", "n3"): [(("met", "m2"), 0.7), (("met", "m3"), 0.4),
                         (("met", "m5"), 1.1)]},
        # a CYCLE: two added reactions sharing two metabolites
        {("rxn", "cy1"): [(("met", "m1"), 1.0), (("met", "m2"), 1.0)],
         ("rxn", "cy2"): [(("met", "m1"), 2.0), (("met", "m2"), 3.0)]},
        # reinforcement: a parallel copy of an edge the host already carries
        {("rxn_reinf", "r0"): [(("met", "m0"), 0.9), (("met", "m1"), 0.9)]},
    ]
    worst_dw = worst_dd = 0.0
    for k, ar2m in enumerate(cases):
        ar2m_v, _, _ = _norm_ar2m(solver, ar2m)
        d_direct = reff_direct(dsolver, ar2m_v, pairs, rb_direct)[0][0]
        d_wood, _, _, _ = reff_summary(solver, ar2m)

        Ga = G.copy()
        for r, mets in ar2m.items():
            for m, c in mets:
                if Ga.has_edge(r, m):
                    Ga[r][m]["w"] = Ga[r][m].get("w", 1.0) + float(c)
                else:
                    Ga.add_edge(r, m, w=float(c))
        for u, v, dd in Ga.edges(data=True):
            dd.setdefault("w", 1.0)
        d_dense = _reff_dense(G, s, t) - _reff_dense(Ga, s, t, wk="w")

        e_dw = abs(d_direct - d_wood)
        e_dd = abs(d_direct - d_dense)
        worst_dw = max(worst_dw, e_dw)
        worst_dd = max(worst_dd, e_dd)
        print(f"case {k}: direct={d_direct:.12f} woodbury={d_wood:.12f} "
              f"dense={d_dense:.12f} | d-w {e_dw:.2e}  d-dense {e_dd:.2e}")
        assert d_direct >= 0.0, "delta_reff < 0 from DirectSolver"

    # every axis at once must equal the axis-at-a-time answer -- this is the
    # property that lets the axis loop die, so it is asserted, not assumed.
    multi = SMWSolver(G, [("met", "m0"), ("met", "m2")],
                      [("met", "m5"), ("met", "m4")],
                      edge_weight_key=None, context=ctx)
    rbm = dsolver.base_reff(multi.pairs)
    for k in range(len(multi.pairs)):
        assert abs(rbm[k] - reff_base(multi, k)) < 1e-9, \
            f"multi-axis base mismatch at pair {k} -- the pair=0 bug"
    print(f"multi-axis: {len(multi.pairs)} pairs agree with per-pair Woodbury")

    assert worst_dw < 1e-9 and worst_dd < 1e-9, \
        f"direct disagrees (vs woodbury {worst_dw:.2e}, vs dense {worst_dd:.2e})"
    print(f"\nPASS -- direct vs woodbury {worst_dw:.2e}, vs dense {worst_dd:.2e}")
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if cmd == "selftest":
        sys.exit(_selftest())
    if cmd == "selftest-direct":
        sys.exit(_selftest_direct())
    print(f"unknown subcommand: {cmd}", file=sys.stderr)
    sys.exit(2)
