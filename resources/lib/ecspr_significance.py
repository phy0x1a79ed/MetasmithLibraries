"""ECSPr significance -- FULL log-normal MIXTURE survival function vs a staged null.

Port of the canonical scadc scorer (`04_reaction_network/reff/05r_significance_mix.py`
and its `_stats.py` dep). `selftest` reproduces that scorer's table to a pinned
tolerance and exits non-zero otherwise -- it is CPU-only and runs in seconds, so it
can gate every commit.

THE SCORER
----------
The p-value answers: *how often does a random insert of the same size match or beat
this Delta?* A GMM is fitted to log(Delta>0) per cell and the upper tail of the WHOLE
mixture is integrated -- every component, at its fitted weight:

    p(x) = sum_k w_k * Phibar((ln x - mu_k) / sigma_k)     for x > 0;  p = 1 otherwise

The GMM is a DENSITY ESTIMATOR here, not a partition tool. Earlier scorers tested
against the "does-nothing" background component alone, discarding the null's higher
random-hit modes -- but those modes are things random draws genuinely do, so excluding
them thins the null and manufactures significance for cells sitting inside them.

SIZE MATCHING
-------------
Observations are size-matched by ORF count to the two flanking null anchors and their
SURVIVAL FUNCTIONS are interpolated:

    p(x) = (1-w) * S_lo(x) + w * S_hi(x)

A convex combination of SFs is itself an SF -- that of the proximity-weighted mixture
null, i.e. a null generated as-if at the fosmid's true size. Interpolating fitted
PARAMETERS instead has no such interpretation; do not "simplify" this into a lerp of
mu/sigma. Nearest-anchor snapping (the retired approach) is a step artifact, since the
null scale grows substantially across the anchor range.

THE TAIL IS HONEST ONLY AS FAR AS THE DRAWS GO
----------------------------------------------
`p_floor` = 1/(n_nonzero+1) marks where counting resolution ends. `p` is left UNFLOORED
so out-of-distribution cells get an estimate rather than a bound; `at_floor` flags rows
whose p rests on the FITTED TAIL rather than on counted draws. Those rows are
extrapolation and must be labelled as such -- they carry no headline. Every row also
carries `p_emp`: pure Phipson & Smyth counting on the SAME raw pool. p<alpha in both =
supported by counting; mixture-only = extrapolated.

DELIBERATELY NOT DONE
---------------------
  * p is NOT scaled by (1-pi0). It is conditional on Delta>0, matching house policy.
    We only report rows with delta_obs>0, so the conditional tail is the correct p.
  * The fit is on the RAW cell draws INCLUDING structural zeros. The fitter derives the
    zero mass itself; feeding it only nonzeros silently zeroes pi0 and changes every p.
  * There is NO DRAW_SIZES constant in this module, by design. Sizes are DERIVED from
    the staged nulls directory. A constant here forked from the null generator's copy
    once already, inside a single directory, with nothing failing. The caller stages a
    curated directory built from an explicit list; this module trusts that directory and
    nothing else. Never point it at a raw cache -- caches accumulate retired sizes.
  * kmax is PINNED (--kmax) rather than left implicit: the fitter's own default is lower
    than the canonical cap, so an implicit cap would make the comparison dishonest.
  * q is computed for schema parity only; ranking is by effect size (delta_obs), which
    is per-row and basis-independent.

Env: numpy + pandas + scikit-learn + scipy (CPU). CLI: selftest | run.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

Q_THRESHOLD = 0.05
ELEMENTS = ["C", "N", "S", "P"]
STYLES = ["A", "B", "D", "E"]
EPS = 1e-12
KMAX = 4            # pinned; the fitter's default is lower -- see module docstring
SEED = 0
LANES = {"reff": ("delta_reff", "reff", "effective resistance"),
         "ieff": ("delta_ieff", "ieff", "effective current (conductance)")}

_NULL_RE = re.compile(r"^(?P<prefix>\w+)_null_canonical_N(?P<n>\d+)\.tsv$")


# =====================================================================
# vendored verbatim from scadc `_stats.py` -- do not "improve" in place
# =====================================================================
# These are byte-faithful copies of the canonical helpers. The parity gate is only
# meaningful if they are. A previous copy of `fit_lognormal_mixture` here had drifted
# to a PRE-GUARD revision (lowest-mu substantial component rather than the guarded
# largest-weight one); it happened not to move the mixture SF -- which reads
# `components`, populated before the guard runs -- but a gate resting on "happened not
# to" is not a gate. If scadc's `_stats.py` changes, re-copy; do not hand-merge.

def bh_q(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def _gmm_icl(gm, Y) -> float:
    """ICL = BIC + 2*EN, EN = -sum_i sum_k tau_ik ln tau_ik (responsibility entropy).

    The entropy term penalizes components whose responsibilities overlap, so a smooth
    heavy tail is NOT split into near-duplicate Gaussians the way plain BIC does at
    large n. For K=1 every tau=1 and EN=0, so ICL == BIC (no penalty change).
    """
    tau = gm.predict_proba(Y)
    t = tau[tau > 0.0]                    # mask zeros so log(0) is never evaluated
    en = -float(np.sum(t * np.log(t)))
    return float(gm.bic(Y) + 2.0 * en)


def fit_lognormal_mixture(null_deltas, kmax: int = 3, min_weight: float = 0.05,
                          min_sigma: float = 1e-3, seed: int = 0) -> dict:
    """Fit a Gaussian mixture to log(Delta>0) and pick the lowest-Delta component.

    Component count K in [1, kmax] chosen by ICL. Background selection is GUARDED:
    the background is the largest-weight component, accepted only if it is (a) a
    majority and (b) the furthest-left / lowest-mu component; otherwise the split is
    untrusted and the background falls back to a single log-normal over ALL nonzero
    draws (bg_all=True), which is conservative.

    NOTE for this module: the mixture SF does not use the background selection at all
    -- it reads `components`, which is populated before the guard runs and is
    well-defined on every path. The guard is retained because this is a verbatim copy.
    """
    out = {"mu_low": np.nan, "sigma_low": np.nan, "w_low": np.nan,
           "k": 0, "icl": np.nan, "bic": np.nan, "pi0": np.nan, "n_nonzero": 0,
           "separated": False, "higher": [], "components": [], "bg_all": False,
           "flag": ""}
    if null_deltas is None or len(null_deltas) == 0:
        out["flag"] = "empty"
        return out
    d = np.asarray(null_deltas, dtype=float)
    out["pi0"] = float((d <= 0.0).mean())
    y = np.log(d[d > 0.0])
    out["n_nonzero"] = int(len(y))

    # too few nonzero draws to fit a mixture: single log-normal MLE
    if len(y) < 20:
        if len(y) == 0:
            out["flag"] = "all-zero"
            return out
        mu1 = float(y.mean())
        sig1 = float(max(y.std(ddof=1) if len(y) > 1 else min_sigma, min_sigma))
        out.update(mu_low=mu1, sigma_low=sig1, w_low=1.0, k=1, separated=False,
                   components=[{"weight": 1.0, "mu": mu1, "sigma": sig1}],
                   flag="too-few")
        return out

    try:
        from sklearn.mixture import GaussianMixture
    except Exception:
        out.update(mu_low=float(y.mean()),
                   sigma_low=float(max(y.std(ddof=1), min_sigma)),
                   w_low=1.0, k=1, flag="no-sklearn")
        return out

    Y = y.reshape(-1, 1)
    best = None  # (gm, icl, k)
    for k in range(1, kmax + 1):
        if k > len(y):
            break
        try:
            gm = GaussianMixture(n_components=k, covariance_type="full",
                                 random_state=seed, max_iter=200, n_init=2,
                                 reg_covar=1e-6)
            gm.fit(Y)
            icl = _gmm_icl(gm, Y)
        except Exception:
            continue
        if best is None or icl < best[1]:
            best = (gm, icl, k)
    if best is None:
        out.update(mu_low=float(y.mean()),
                   sigma_low=float(max(y.std(ddof=1), min_sigma)),
                   w_low=1.0, k=1, flag="fit-failed")
        return out

    gm, icl, k = best
    bic = float(gm.bic(Y))
    mus = gm.means_.ravel()
    sigmas = np.sqrt(gm.covariances_.ravel())
    ws = gm.weights_.ravel()

    order = np.argsort(mus)
    out["components"] = [
        {"weight": float(ws[i]), "mu": float(mus[i]),
         "sigma": float(max(sigmas[i], min_sigma))} for i in order]
    out.update(k=int(k), icl=float(icl), bic=float(bic))

    L = int(np.argmax(ws))
    majority = ws[L] >= 0.5
    furthest_left = (L == int(order[0]))
    if not (majority and furthest_left):
        mu_all = float(y.mean())
        sig_all = float(max(y.std(ddof=1) if len(y) > 1 else min_sigma, min_sigma))
        out.update(mu_low=mu_all, sigma_low=sig_all, w_low=1.0, bg_all=True,
                   separated=False, higher=[],
                   flag=("k1" if k == 1 else "guard-fallback"))
        return out

    out.update(mu_low=float(mus[L]),
               sigma_low=float(max(sigmas[L], min_sigma)),
               w_low=float(ws[L]), bg_all=False, separated=bool(k > 1))
    out["higher"] = [
        {"weight": float(ws[i]), "mu": float(mus[i]),
         "median": float(np.exp(mus[i]))}
        for i in order if mus[i] > mus[L]]
    if k == 1:
        out["flag"] = "k1"
    return out


# =====================================================================
# scorers
# =====================================================================

def mixture_scorer(d: np.ndarray, kmax: int = KMAX, seed: int = SEED):
    """Upper-tail area of the FULL fitted log-normal mixture over the raw draws.

    `d` MUST be the raw cell draws including structural zeros -- the fitter derives
    pi0 from (d <= 0) itself and fits log(d[d>0]).

    Returns (p_fn, flag, n_nonzero, p_floor, fit). The SF is left unfloored so OOD
    cells get an estimate; p_floor only records where counting resolution ends.
    """
    fit = fit_lognormal_mixture(np.asarray(d, dtype=float), kmax=kmax, seed=seed)
    comps = fit.get("components") or []
    if not comps:
        # no-sklearn / fit-failed: the fitter still returns a finite single log-normal
        # in mu_low/sigma_low. Fall back to it rather than dropping the cell.
        mu, sg = fit.get("mu_low", np.nan), fit.get("sigma_low", np.nan)
        if np.isfinite(mu) and np.isfinite(sg) and sg > 0:
            comps = [{"weight": 1.0, "mu": float(mu), "sigma": float(sg)}]
    n = int(fit.get("n_nonzero", 0))
    if not comps or n == 0:
        return (lambda x: 1.0), (fit.get("flag") or "degenerate"), 0, np.nan, fit

    w = np.array([c["weight"] for c in comps], dtype=float)
    w = w / w.sum()                                   # GMM weights already sum to 1
    mu = np.array([c["mu"] for c in comps], dtype=float)
    sg = np.maximum(np.array([c["sigma"] for c in comps], dtype=float), 1e-12)
    floor = 1.0 / (n + 1.0)

    def p(x):
        if x is None or x <= 0:
            return 1.0
        return float(np.clip(float(w @ norm.sf((np.log(x) - mu) / sg)), 1e-300, 1.0))
    return p, fit.get("flag", ""), n, floor, fit


def empirical_scorer(d: np.ndarray):
    """Phipson & Smyth exceedance on the SAME raw pool -- the counting reference.

    p(x) = (1 + #{nz >= x}) / (1 + n_nz). Fed the full pool (not a background split),
    so `p_emp` isolates fitted-tail extrapolation from counted evidence within one row.
    """
    nz = np.asarray(d, dtype=float)
    nz = np.sort(nz[nz > 0])
    n = len(nz)
    if n == 0:
        return lambda x: 1.0

    def p(x):
        if x is None or x <= 0:
            return 1.0
        return float((1 + (n - int(np.searchsorted(nz, x, side="left")))) / (n + 1))
    return p


# =====================================================================
# null loading -- sizes DERIVED from the staged directory, never a constant
# =====================================================================

def discover_draw_sizes(nulls_dir, prefix: str) -> list[int]:
    """Sorted draw sizes present in the STAGED nulls directory.

    This is the module's only source of draw sizes. The directory must be the curated
    one the driver builds from an explicit list -- a raw cache also holds retired sizes
    and, in at least one case, an un-suffixed leftover from an overwrite incident.
    """
    sizes = sorted(
        int(m.group("n"))
        for f in Path(nulls_dir).iterdir()
        if (m := _NULL_RE.match(f.name)) and m.group("prefix") == prefix
    )
    if not sizes:
        raise SystemExit(
            f"no {prefix}_null_canonical_N*.tsv under {nulls_dir}; the staged nulls "
            f"directory is the only source of draw sizes for this module"
        )
    return sizes


def bracket(orf_count: int, draw_sizes: list[int]):
    """Two anchors flanking orf_count + interpolation weight w in [0,1].

    Clamped outside the anchor range. See module docstring: the two anchors' SURVIVAL
    FUNCTIONS are interpolated, never their parameters.
    """
    lo, hi = draw_sizes[0], draw_sizes[-1]
    n = min(max(int(orf_count), lo), hi)
    for a, b in zip(draw_sizes[:-1], draw_sizes[1:]):
        if a <= n <= b:
            return a, b, (0.0 if b == a else (n - a) / (b - a))
    return hi, hi, 0.0


# =====================================================================
# ORF-count sources
# =====================================================================

def fosmid_orf_counts_from_faa(path) -> dict:
    """Count ORFs per contig from a protein FASTA (>contig_<geneN> ...).

    contig id = ORF id minus the trailing _<geneN>, via rsplit. Split-contig ids like
    `C00310.A` survive -- a `\\w+`-based regex would DROP them silently, because `.` is
    not a word character, and the fosmid would vanish from the output with no error.
    """
    counts: dict = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                counts[line[1:].split()[0].rsplit("_", 1)[0]] += 1
    return dict(counts)


def fosmid_orf_counts_from_evidence(ev_path, source: str = "fosmid") -> dict:
    """Distinct ORFs per contig from the compiled evidence parquet."""
    ev = pd.read_parquet(ev_path, columns=["source", "orf"])
    ev = ev[ev["source"] == source].copy()
    ev["contig"] = ev["orf"].str.rsplit("_", n=1).str[0]
    return ev.groupby("contig")["orf"].nunique().astype(int).to_dict()


# =====================================================================
# scoring
# =====================================================================

def load_nulls(nulls_dir, metric: str, prefix: str, draw_sizes: list[int]) -> dict:
    nulls_dir = Path(nulls_dir)
    nul = {}
    for N in draw_sizes:
        p = nulls_dir / f"{prefix}_null_canonical_N{N}.tsv"
        if not p.exists():
            raise SystemExit(f"missing staged null {p}")
        df = pd.read_csv(p, sep="\t", usecols=["element", "axis_id", "null", metric])
        df.loc[df[metric] < EPS, metric] = 0.0
        nul[N] = df[df["element"].isin(ELEMENTS)]
    return nul


def score(nul: dict, obs: pd.DataFrame, metric: str, draw_sizes: list[int],
          kmax: int = KMAX, seed: int = SEED, verbose: bool = True) -> pd.DataFrame:
    cell_p, cell_e, cell_n, cell_floor, cell_k = {}, {}, {}, {}, {}
    for N in draw_sizes:
        flags: dict = defaultdict(int)
        for (el, ax, st), g in nul[N].groupby(["element", "axis_id", "null"],
                                              sort=False):
            if st not in STYLES:
                continue
            v = g[metric].values.astype(float)      # RAW pool: zeros included
            fn, flag, n, floor, fit = mixture_scorer(v, kmax=kmax, seed=seed)
            cell_p[(el, ax, st, N)] = fn
            cell_e[(el, ax, st, N)] = empirical_scorer(v)
            cell_n[(el, ax, st, N)] = n
            cell_floor[(el, ax, st, N)] = floor
            cell_k[(el, ax, st, N)] = int(fit.get("k", 0))
            flags[flag or "ok"] += 1
        if verbose:
            ks = [cell_k[k] for k in cell_k if k[3] == N]
            print(f"  N={N}: {len(ks)} cells fitted | ICL k: "
                  + " ".join(f"k={i}:{ks.count(i)}" for i in sorted(set(ks)))
                  + " | fitter flags: "
                  + " ".join(f"{k}:{v}" for k, v in sorted(flags.items())), flush=True)

    recs = []
    for r in obs.itertuples():
        Nlo, Nhi, w = bracket(r.n_orfs, draw_sizes)
        for st in STYLES:
            klo = (r.element, r.axis_id, st, Nlo)
            khi = (r.element, r.axis_id, st, Nhi)
            fnlo, fnhi = cell_p.get(klo), cell_p.get(khi)
            if fnlo is None or fnhi is None:
                continue
            x = float(getattr(r, metric))
            # Interpolate the two anchors' SF VALUES, never their (mu, sigma).
            p = (1.0 - w) * fnlo(x) + w * fnhi(x)
            pe = (1.0 - w) * cell_e[klo](x) + w * cell_e[khi](x)
            pfloor = (1.0 - w) * cell_floor[klo] + w * cell_floor[khi]
            recs.append(dict(fosmid=r.fosmid, element=r.element,
                             axis_id=r.axis_id, n_orfs=int(r.n_orfs),
                             N_lo=Nlo, N_hi=Nhi, w=round(float(w), 3), null=st,
                             delta_obs=x, p=p, p_emp=pe,
                             k_lo=cell_k[klo], k_hi=cell_k[khi],
                             n_bg=int(min(cell_n[klo], cell_n[khi])),
                             p_floor=pfloor,
                             at_floor=bool(x > 0 and p < pfloor)))
    sig = pd.DataFrame(recs)
    if sig.empty:
        sig["q"] = []
        sig["survives"] = []
        return sig
    # q for schema parity only -- ranking is by effect size (delta_obs).
    sig["q"] = 1.0
    for (_el, _st), g in sig.groupby(["element", "null"], sort=False):
        tested = g[(g["delta_obs"] > 0) & g["p"].notna()]
        if len(tested):
            sig.loc[tested.index, "q"] = bh_q(tested["p"].values)
    sig["survives"] = sig["q"] < Q_THRESHOLD
    return sig


def _load_obs(report_path, metric: str) -> pd.DataFrame:
    obs = pd.read_csv(report_path, sep="\t",
                      usecols=["element", "axis_id", "fosmid", "null", metric])
    obs = obs[(obs["null"] == "obs") & obs["element"].isin(ELEMENTS)].copy()
    obs.loc[obs[metric] < EPS, metric] = 0.0
    return obs


def _prepare(report, nulls_dir, metric, prefix, orf_counts, kmax, seed,
             verbose=True):
    draw_sizes = discover_draw_sizes(nulls_dir, prefix)
    if verbose:
        print(f"[sig:{prefix}] draw sizes from staged nulls: {draw_sizes}", flush=True)
    nul = load_nulls(nulls_dir, metric, prefix, draw_sizes)
    obs = _load_obs(report, metric)
    obs["n_orfs"] = obs["fosmid"].astype(str).map(orf_counts)
    missing = obs[obs["n_orfs"].isna()]["fosmid"].nunique()
    obs = obs[obs["n_orfs"].notna()].copy()
    obs["n_orfs"] = obs["n_orfs"].astype(int)
    if verbose:
        print(f"[sig:{prefix}] {obs['fosmid'].nunique()} fosmids x "
              f"{obs['axis_id'].nunique()} axes ({len(obs)} obs rows), kmax={kmax}"
              + (f" | {missing} fosmid(s) had no ORF count and were dropped"
                 if missing else ""), flush=True)
    return nul, obs, draw_sizes


def run(args):
    metric, prefix, _lane_name = LANES[args.lane]
    if args.faa:
        fc = fosmid_orf_counts_from_faa(args.faa)
    elif args.evidence:
        fc = fosmid_orf_counts_from_evidence(args.evidence)
    else:
        raise SystemExit("need --faa or --evidence for ORF counts")

    nul, obs, draw_sizes = _prepare(args.report, args.nulls_dir, metric, prefix,
                                    fc, args.kmax, args.seed)
    sig = score(nul, obs, metric, draw_sizes, kmax=args.kmax, seed=args.seed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(args.out, sep="\t", index=False)
    if not sig.empty:
        pos = sig[sig.delta_obs > 0]
        nfloor = int(pos["at_floor"].sum())
        print(f"[sig:{prefix}] {int((pos.p < 0.05).sum())} cells at raw p<0.05 "
              f"({int(((pos.p < 0.05) & (pos.p_emp < 0.05)).sum())} supported by "
              f"counting, rest extrapolated) | {nfloor}/{len(pos)} positive-delta "
              f"rows rest on the fitted tail", flush=True)
    print(f"[sig:{prefix}] wrote {args.out}", flush=True)


# =====================================================================
# selftest: reproduce the canonical scadc table
# =====================================================================

def selftest(args):
    """Join against the canonical scadc scorer's table and fail past the tolerance.

    Compares `p` and `p_emp` numerically and `at_floor` exactly. Also asserts the
    split-contig ids survive the ORF counter -- their absence is the `\\w`-regex bug,
    and it is silent, so it has to be checked positively rather than noticed.
    """
    ref_path = Path(args.reference)
    nulls_dir = Path(args.nulls_dir)
    report_dir = Path(args.report_dir)
    fc = fosmid_orf_counts_from_faa(args.faa)

    metric, prefix, _ = LANES[args.lane]
    nul, obs, draw_sizes = _prepare(report_dir / f"{prefix}_axes_report.tsv",
                                    nulls_dir, metric, prefix, fc,
                                    args.kmax, args.seed)
    sig = score(nul, obs, metric, draw_sizes, kmax=args.kmax, seed=args.seed)
    ref = pd.read_csv(ref_path, sep="\t")

    # Refuse a stale reference by SHAPE before comparing, so pointing the gate at a
    # retired sibling table fails with a sentence instead of a KeyError three frames
    # deep in pandas. `matched_N` belonged to the retired nearest-size scorer;
    # `p_emp` did not exist until the mixture scorer emitted it per row.
    if "matched_N" in ref.columns:
        raise SystemExit(
            f"[selftest:{prefix}] reference {ref_path.name} carries `matched_N` -- that "
            f"is the retired nearest-size scorer's schema, not the full-mixture SF's. "
            f"This gate reproduces the canonical table; point --reference at it.")
    ref_missing = [c for c in ("p", "p_emp", "at_floor") if c not in ref.columns]
    if ref_missing:
        raise SystemExit(
            f"[selftest:{prefix}] reference {ref_path.name} lacks {ref_missing}; it is "
            f"not the canonical scorer's output. A retired sibling table has no `p_emp` "
            f"column because the empirical null was never a competing table -- it is a "
            f"column on every canonical row.")

    keys = ["fosmid", "element", "axis_id", "null"]
    m = sig.merge(ref, on=keys, suffixes=("_mine", "_ref"))
    if not len(m):
        raise SystemExit(f"[selftest:{prefix}] FAILED: no rows joined against {ref_path}")

    dp = float(np.abs(m["p_mine"] - m["p_ref"]).max())
    dpe = float(np.abs(m["p_emp_mine"] - m["p_emp_ref"]).max())
    floor_disagree = int((m["at_floor_mine"] != m["at_floor_ref"]).sum())

    print(f"[selftest:{prefix}] rows mine={len(sig)} ref={len(ref)} matched={len(m)}")
    print(f"[selftest:{prefix}] max|dp|={dp:.3e}  max|dp_emp|={dpe:.3e}  "
          f"at_floor disagreements={floor_disagree}")

    # positive check: the split contigs must be present, not merely un-errored
    mine_fos = set(sig["fosmid"].astype(str))
    missing = [c for c in args.split_contigs if c not in mine_fos]
    if missing:
        raise SystemExit(
            f"[selftest:{prefix}] FAILED: split contigs absent from output: {missing}. "
            f"A `\\w`-based ORF-id regex drops these silently.")
    print(f"[selftest:{prefix}] split contigs present: {list(args.split_contigs)}")

    unmatched = len(ref) - len(m)
    if unmatched:
        print(f"[selftest:{prefix}] WARNING: {unmatched} reference row(s) did not join")

    worst = max(dp, dpe)
    if worst > args.tol or floor_disagree or unmatched:
        raise SystemExit(
            f"[selftest:{prefix}] PARITY FAILED: worst|d|={worst:.3e} "
            f"(tol {args.tol:.1e}), at_floor disagreements={floor_disagree}, "
            f"unjoined ref rows={unmatched}")
    print(f"[selftest:{prefix}] PASS (tol {args.tol:.1e})")


def parse_args():
    ap = argparse.ArgumentParser(
        description="ECSPr significance (full-mixture SF vs a staged null)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    for name, help_ in (("run", "score an observed report vs the staged null"),
                        ("selftest", "reproduce the canonical scadc table")):
        p = sub.add_parser(name, help=help_)
        p.add_argument("--lane", choices=list(LANES), required=True)
        p.add_argument("--nulls-dir", required=True,
                       help="CURATED staged nulls dir; draw sizes are derived from it")
        p.add_argument("--kmax", type=int, default=KMAX,
                       help="ICL component-count cap (pinned; see module docstring)")
        p.add_argument("--seed", type=int, default=SEED)

    run_p, self_p = sub.choices["run"], sub.choices["selftest"]

    run_p.add_argument("--report", required=True, help="ecspr_network solve output tsv")
    run_p.add_argument("--faa", default=None, help="proteins (>contig_N) for ORF counts")
    run_p.add_argument("--evidence", default=None, help="evidence parquet fallback")
    run_p.add_argument("--out", required=True)
    run_p.set_defaults(func=run)

    self_p.add_argument("--reference", required=True,
                        help="canonical scadc significance table to join against")
    self_p.add_argument("--report-dir", required=True,
                        help="dir holding {lane}_axes_report.tsv")
    self_p.add_argument("--faa", required=True)
    self_p.add_argument("--tol", type=float, required=True,
                        help="pinned by the caller; do not relax to pass")
    self_p.add_argument("--split-contigs", nargs="*", default=[],
                        help="ids that must survive the ORF counter")
    self_p.set_defaults(func=selftest)
    return ap.parse_args()


if __name__ == "__main__":
    a = parse_args()
    a.func(a)
