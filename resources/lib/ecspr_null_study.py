"""ECSPr null-family STUDY scorer -- try a different null model, safely.

Null GENERATION (how draws are made: style, K, N, seed) and null SCORING (how a
p-value comes off those draws) are orthogonal. This module is the scoring half: it
takes the same staged draws the canonical scorer uses and applies a DIFFERENT
family to them, so a new null model is a question you can ask without forking the
pipeline that answers it.

WHY THIS IS A SEPARATE OUTPUT AND NOT A FLAG
--------------------------------------------
This writes `ecspr::{lane}_significance_study`, never `ecspr::{lane}_significance`.
The canonical significance type keeps exactly ONE producer -- the full-mixture SF --
and no switch can change what it means.

That is deliberate and it is the lesson of the incumbent. When the scorer became a
swappable flag, three figure generators drifted onto three different null tables
across two days, and one shipped a title asserting the opposite of its own
annotations. A default-able family flag is the mechanism that produced that. So an
alternative family gets a separate, clearly-labelled table that no figure reads by
accident, and promoting one is a deliberate edit to the canonical scorer, reviewed
and re-gated -- not a parameter someone passes.

The family is a STAGED, CONTENT-HASHED input, so two families produce different
cache keys and cannot collide on one entry. That is the same property that makes
the compute profile safe, applied to the science knob it was invented for.

FAMILIES
--------
  mixture    the canonical full-mixture SF. Present as the REFERENCE, so a study
             run can show its family against the incumbent on one axis, in one
             table, from one pool.
  empirical  Phipson & Smyth counting on the raw pool. Exactly invariant under any
             monotone reparametrisation (counting is), hence the honest ceiling --
             but it cannot resolve past 1/(n+1).
  gpd        zero-inflated peaks-over-threshold with a generalized-Pareto tail.
             A genuinely different family: it extrapolates past the draws where
             counting stops. A Vuong test favours a power-law tail over log-normal
             in nearly every cell, and the two disagree by orders of magnitude far
             out -- which is the open question a study table exists to explore.

Adding a family = one entry in FAMILIES. Nothing else changes.

Env: numpy + pandas + scikit-learn + scipy (CPU). CLI: run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from ecspr_significance import (  # noqa: E402
    ELEMENTS, EPS, KMAX, LANES, SEED, STYLES,
    _load_obs, bh_q, bracket, discover_draw_sizes, empirical_scorer,
    fit_lognormal_mixture, fosmid_orf_counts_from_faa, load_nulls, mixture_scorer,
)

Q_THRESHOLD = 0.05


# =====================================================================
# families
# =====================================================================

def _family_mixture(d, kmax, seed):
    fn, _flag, n, floor, _fit = mixture_scorer(d, kmax=kmax, seed=seed)
    return fn, n, floor


def _family_empirical(d, kmax, seed):
    nz = np.asarray(d, dtype=float)
    n = int((nz > 0).sum())
    return empirical_scorer(d), n, (1.0 / (n + 1.0) if n else np.nan)


def _family_gpd(d, kmax, seed, tail_q: float = 0.80):
    """Zero-inflated POT-GPD, conditional on Delta>0.

    Conditional to match the other families -- the pi0 factor is divided back out.
    Without that the comparison is not like-for-like and the family looks better
    than it is by exactly 1/(1-pi0).
    """
    from scipy.stats import genpareto
    d = np.asarray(d, dtype=float)
    nz = np.sort(d[d > 0])
    n = len(nz)
    if n == 0:
        return (lambda x: 1.0), 0, np.nan
    floor = 1.0 / (n + 1.0)
    if n < 20:
        return empirical_scorer(d), n, floor

    u = float(np.quantile(nz, tail_q))
    exc = nz[nz > u] - u
    frac_above_u = len(exc) / float(n)
    fit = None
    if len(exc) >= 8:
        try:
            xi, _loc, sigma = genpareto.fit(exc, floc=0.0)
            if np.isfinite(xi) and np.isfinite(sigma) and sigma > 0:
                fit = (xi, sigma)
        except Exception:
            fit = None

    def p(x):
        if x is None or x <= 0:
            return 1.0
        if x <= u:
            return float((np.sum(nz >= x) + 1.0) / (n + 1.0))
        if fit is None:
            return float(max(frac_above_u * np.exp(-(x - u) / max(exc.mean(), 1e-12)),
                             1e-300))
        xi, sigma = fit
        sf = float(genpareto.sf(x - u, xi, loc=0.0, scale=sigma))
        return float(np.clip(frac_above_u * max(sf, 0.0), 1e-300, 1.0))
    return p, n, floor


FAMILIES = {
    "mixture": _family_mixture,
    "empirical": _family_empirical,
    "gpd": _family_gpd,
}


def read_family_spec(path) -> str:
    """`family: <name>` -- a staged, content-hashed input, not a CLI default.

    Fails loudly on an unknown family. A study that silently fell back to the
    canonical family would be indistinguishable from a study that confirmed it.
    """
    fam = None
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        k, _, v = line.partition(":")
        if k.strip() == "family":
            fam = v.strip()
    if fam not in FAMILIES:
        raise SystemExit(
            f"null_model_spec: family {fam!r} not in {sorted(FAMILIES)}")
    return fam


# =====================================================================
# scoring
# =====================================================================

def score_family(nul, obs, metric, draw_sizes, family: str, kmax=KMAX, seed=SEED):
    fn_family = FAMILIES[family]
    cell_p, cell_n, cell_floor = {}, {}, {}
    for N in draw_sizes:
        for (el, ax, st), g in nul[N].groupby(["element", "axis_id", "null"],
                                              sort=False):
            if st not in STYLES:
                continue
            v = g[metric].values.astype(float)      # RAW pool: zeros included
            fn, n, floor = fn_family(v, kmax, seed)
            cell_p[(el, ax, st, N)] = fn
            cell_n[(el, ax, st, N)] = n
            cell_floor[(el, ax, st, N)] = floor

    recs = []
    for r in obs.itertuples():
        Nlo, Nhi, w = bracket(r.n_orfs, draw_sizes)
        for st in STYLES:
            klo, khi = (r.element, r.axis_id, st, Nlo), (r.element, r.axis_id, st, Nhi)
            fnlo, fnhi = cell_p.get(klo), cell_p.get(khi)
            if fnlo is None or fnhi is None:
                continue
            x = float(getattr(r, metric))
            # Same SF interpolation as the canonical scorer: a convex combination of
            # SFs is an SF. Families differ in the SF, never in how sizes are matched
            # -- otherwise a study measures two changes at once and attributes both
            # to the family.
            p = (1.0 - w) * fnlo(x) + w * fnhi(x)
            pfloor = (1.0 - w) * cell_floor[klo] + w * cell_floor[khi]
            recs.append(dict(fosmid=r.fosmid, element=r.element, axis_id=r.axis_id,
                             n_orfs=int(r.n_orfs), N_lo=Nlo, N_hi=Nhi,
                             w=round(float(w), 3), null=st, delta_obs=x,
                             family=family, p=p,
                             n_bg=int(min(cell_n[klo], cell_n[khi])),
                             p_floor=pfloor,
                             at_floor=bool(x > 0 and p < pfloor)))
    sig = pd.DataFrame(recs)
    if sig.empty:
        sig["q"] = []
        return sig
    sig["q"] = 1.0
    for (_el, _st), g in sig.groupby(["element", "null"], sort=False):
        tested = g[(g["delta_obs"] > 0) & g["p"].notna()]
        if len(tested):
            sig.loc[tested.index, "q"] = bh_q(tested["p"].values)
    return sig


def run(args):
    metric, prefix, _ = LANES[args.lane]
    family = read_family_spec(args.family_spec) if args.family_spec else args.family
    if family not in FAMILIES:
        raise SystemExit(f"family {family!r} not in {sorted(FAMILIES)}")

    fc = fosmid_orf_counts_from_faa(args.faa)
    draw_sizes = discover_draw_sizes(args.nulls_dir, prefix)
    print(f"[study:{prefix}] family={family} | draw sizes from staged nulls: "
          f"{draw_sizes}", flush=True)
    nul = load_nulls(args.nulls_dir, metric, prefix, draw_sizes)
    obs = _load_obs(args.report, metric)
    obs["n_orfs"] = obs["fosmid"].astype(str).map(fc)
    obs = obs[obs["n_orfs"].notna()].copy()
    obs["n_orfs"] = obs["n_orfs"].astype(int)
    print(f"[study:{prefix}] {obs['fosmid'].nunique()} fosmids x "
          f"{obs['axis_id'].nunique()} axes", flush=True)

    sig = score_family(nul, obs, metric, draw_sizes, family,
                       kmax=args.kmax, seed=args.seed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(args.out, sep="\t", index=False)
    pos = sig[sig.delta_obs > 0]
    print(f"[study:{prefix}] family={family}: {int((pos.p < 0.05).sum())} cells at "
          f"raw p<0.05 | {int(pos.at_floor.sum())}/{len(pos)} on the fitted tail",
          flush=True)
    print(f"[study:{prefix}] wrote {args.out} -- a STUDY table. The canonical "
          f"significance type has one producer and is not this.", flush=True)


def parse_args():
    ap = argparse.ArgumentParser(description="ECSPr null-family study scorer")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--lane", choices=list(LANES), required=True)
    r.add_argument("--report", required=True)
    r.add_argument("--nulls-dir", required=True)
    r.add_argument("--faa", required=True)
    r.add_argument("--family", default=None, choices=sorted(FAMILIES))
    r.add_argument("--family-spec", default=None,
                   help="staged null_model_spec (preferred: it is content-hashed)")
    r.add_argument("--kmax", type=int, default=KMAX)
    r.add_argument("--seed", type=int, default=SEED)
    r.add_argument("--out", required=True)
    r.set_defaults(func=run)
    return ap.parse_args()


if __name__ == "__main__":
    a = parse_args()
    a.func(a)
