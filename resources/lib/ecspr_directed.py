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
                  floor=DIODE_BACKWARD_FLOOR, return_iters=False):
    """Two-terminal effective conductance of a rectified network by semismooth Newton.

    Parameters
    ----------
    B : (m, n) signed incidence (CSR), from :func:`build_incidence`.
    gp, gm : (m,) forward / backward conductances per edge (``gm`` is clipped up to
        ``floor * gp`` to keep the grounded system nonsingular).
    s, t : source / sink node indices (unit current injected s -> t).
    g : grounded node index (its potential row/col is dropped; default 0).

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
    phi = np.zeros(n)

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
    return (ceff, nit) if return_iters else ceff


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


def _selftest():
    print("=" * 68)
    print("ecspr_directed self-tests")
    print("=" * 68)
    _selftest_directed_referent()
    _selftest_directed_parity()
    _selftest_directed_behaviour()
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(_selftest())
