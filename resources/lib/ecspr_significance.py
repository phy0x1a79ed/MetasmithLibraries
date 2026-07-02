"""ECSPr significance -- calibrated empirical-CCDF exceedance vs the REUSED null.

Ports scadc `04_reaction_network/reff/05r_significance_emp.py` (+ its `_stats.py`
and `altfits/_modes.py` deps) onto fabfos. Scores each observed per-fosmid
dR_eff / dI_eff (from `ecspr_network` solve) against the **frozen scadc null**
draws -- which are REUSED verbatim, never recomputed (per the fabfos scope: only
the metagenome null may be reused). See ecspr_network / ecspr_solver.

Scoring spine (byte-identical to scadc):
  * size-match each fosmid to the nearest null size N in {21,34,56} by ORF count;
  * per (element, axis, null-style) cell, split the null's nonzero Delta into its
    does-nothing BACKGROUND mode (GMM lowest-mu component, hard responsibility
    assignment) and score the observation with the conservative finite-null
    exceedance p(x) = (1 + #{bg_nz >= x}) / (1 + n)  [Phipson & Smyth 2010],
    floored at 1/(n+1) (~1e-3 for a 1000-draw null -- finer p needs more draws,
    not a parametric tail, all of which are miscalibrated here);
  * BH-q within (element, null-style); survives = q < Q_THRESHOLD.

Observations beyond all null draws tie at the floor and are ranked by EFFECT SIZE
(delta_obs), per the standing protocol. At 1000 reused draws style-D C/P may show
0 BH-survivors -- report effect-size ranks, not survivor counts.

The frozen nulls (reused) live as `{reff,ieff}_null_canonical_N{21,34,56}.tsv`
with cols element, axis_id, null(style), <metric>. The observed report is
`ecspr_network` solve output (null == "obs"). ORF counts come from the FRESH
prodigal proteins on the 132 inserts (--faa) or the compiled evidence (--evidence).

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

Q_THRESHOLD = 0.05
ELEMENTS = ["C", "N", "S", "P"]
DRAW_SIZES = [21, 34, 56]
STYLES = ["A", "B", "D", "E"]
EPS = 1e-12
LANES = {"reff": ("delta_reff", "reff", "effective resistance"),
         "ieff": ("delta_ieff", "ieff", "effective current (conductance)")}


# --- BH q (verbatim from scadc _stats.bh_q) ---------------------------------

def bh_q(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    pvals = np.asarray(pvals, dtype=float)
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


# --- multimodal log-normal mixture (verbatim from scadc _stats) -------------

def _gmm_icl(gm, Y) -> float:
    """ICL = BIC + 2*EN, EN = -sum tau ln tau (responsibility entropy)."""
    tau = gm.predict_proba(Y)
    t = tau[tau > 0.0]
    en = -float(np.sum(t * np.log(t)))
    return float(gm.bic(Y) + 2.0 * en)


def fit_lognormal_mixture(null_deltas, kmax: int = 3, min_weight: float = 0.05,
                          min_sigma: float = 1e-3, seed: int = 0) -> dict:
    """Fit a Gaussian mixture to log(Delta>0), pick lowest-mode as background."""
    out = {"mu_low": np.nan, "sigma_low": np.nan, "w_low": np.nan,
           "k": 0, "icl": np.nan, "bic": np.nan, "pi0": np.nan, "n_nonzero": 0,
           "separated": False, "higher": [], "components": [], "flag": ""}
    if null_deltas is None or len(null_deltas) == 0:
        out["flag"] = "empty"
        return out
    d = np.asarray(null_deltas, dtype=float)
    out["pi0"] = float((d <= 0.0).mean())
    y = np.log(d[d > 0.0])
    out["n_nonzero"] = int(len(y))

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
    best = None
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
    substantial = [i for i in order if ws[i] >= min_weight]
    low = substantial[0] if substantial else int(order[0])

    out.update(mu_low=float(mus[low]),
               sigma_low=float(max(sigmas[low], min_sigma)),
               w_low=float(ws[low]), k=int(k), icl=float(icl), bic=float(bic),
               separated=bool(k > 1))
    out["higher"] = [
        {"weight": float(ws[i]), "mu": float(mus[i]),
         "median": float(np.exp(mus[i]))}
        for i in order if mus[i] > mus[low]]
    out["components"] = [
        {"weight": float(ws[i]), "mu": float(mus[i]),
         "sigma": float(max(sigmas[i], min_sigma))} for i in order]
    if k == 1:
        out["flag"] = "k1"
    return out


# --- background-mode extraction (verbatim from scadc altfits/_modes) --------

def background_draws(d, kmax: int = 4, seed: int = 0):
    """Split a null cell's nonzero Delta into the background mode by GMM assign."""
    from scipy.stats import norm
    d = np.asarray(d, dtype=float)
    nz = d[d > 0]
    fit = fit_lognormal_mixture(d, kmax=kmax, seed=seed)
    comps = fit.get("components", [])
    info = {"k": fit["k"], "w_low": fit["w_low"], "mu_low": fit["mu_low"],
            "separated": fit["separated"], "flag": fit["flag"],
            "n_nonzero": int(len(nz)),
            "frac_bg": (1.0 if len(nz) else 0.0)}

    if len(nz) == 0 or len(comps) <= 1:
        return nz, info

    y = np.log(nz)
    logresp = np.column_stack([
        np.log(max(c["weight"], 1e-300)) + norm.logpdf(y, c["mu"], c["sigma"])
        for c in comps])
    assign = logresp.argmax(axis=1)
    mus = np.array([c["mu"] for c in comps])
    low_comps = np.where(mus <= fit["mu_low"] + 1e-9)[0]
    bg = nz[np.isin(assign, low_comps)]
    info["frac_bg"] = float(len(bg) / len(nz)) if len(nz) else 0.0
    return bg, info


# --- empirical exceedance scorer (verbatim from scadc 05r) ------------------

def empirical_scorer(d: np.ndarray):
    """Conservative empirical exceedance on nonzero background draws.

    p(x) = (1 + #{bg_nz >= x}) / (1 + n), floored at 1/(n+1)."""
    nz = np.asarray(d, dtype=float)
    nz = nz[nz > 0]
    s = np.sort(nz)
    n = len(s)
    if n == 0:
        return (lambda x: 1.0), "no-nonzero", 0
    floor = 1.0 / (n + 1.0)

    def p(x):
        if x is None or x <= 0:
            return 1.0
        ge = n - int(np.searchsorted(s, x, side="left"))
        return float((1 + ge) / (n + 1))
    return p, "", n, floor


def nearest_n(orf_count: int) -> int:
    return int(min(DRAW_SIZES, key=lambda N: abs(N - orf_count)))


# --- ORF-count sources (fabfos: from the FRESH prodigal proteins) -----------

def fosmid_orf_counts_from_faa(path) -> dict:
    """Count ORFs per contig from prodigal amino-acid FASTA (>contig_N ...)."""
    counts: dict = defaultdict(int)
    pat = re.compile(r"^>(\w+?)_\d+")
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                m = pat.match(line)
                if m:
                    counts[m.group(1)] += 1
    return dict(counts)


def fosmid_orf_counts_from_evidence(ev_path, source: str = "fosmid") -> dict:
    """Distinct ORFs per contig from the compiled evidence parquet."""
    ev = pd.read_parquet(ev_path, columns=["source", "orf"])
    ev = ev[ev["source"] == source].copy()
    ev["contig"] = ev["orf"].str.rsplit("_", n=1).str[0]
    return ev.groupby("contig")["orf"].nunique().astype(int).to_dict()


# --- null loading + scoring (port of scadc 05r) -----------------------------

def load_nulls(nulls_dir, metric: str, prefix: str) -> dict:
    nulls_dir = Path(nulls_dir)
    nul = {}
    for N in DRAW_SIZES:
        p = nulls_dir / f"{prefix}_null_canonical_N{N}.tsv"
        if not p.exists():
            raise SystemExit(f"missing frozen null {p}")
        df = pd.read_csv(p, sep="\t", usecols=["element", "axis_id", "null", metric])
        df.loc[df[metric] < EPS, metric] = 0.0
        nul[N] = df[df["element"].isin(ELEMENTS)]
    return nul


def score(nul: dict, obs: pd.DataFrame, metric: str) -> pd.DataFrame:
    cell_p, cell_n, cell_floor = {}, {}, {}
    for N in DRAW_SIZES:
        for (el, ax, st), g in nul[N].groupby(["element", "axis_id", "null"],
                                              sort=False):
            if st not in STYLES:
                continue
            bg, _info = background_draws(g[metric].values)
            fn, _flag, n, *rest = empirical_scorer(bg)
            cell_p[(el, ax, st, N)] = fn
            cell_n[(el, ax, st, N)] = n
            cell_floor[(el, ax, st, N)] = rest[0] if rest else np.nan

    recs = []
    for r in obs.itertuples():
        for st in STYLES:
            key = (r.element, r.axis_id, st, r.matched_N)
            fn = cell_p.get(key)
            if fn is None:
                continue
            x = float(getattr(r, metric))
            recs.append(dict(fosmid=r.fosmid, element=r.element,
                             axis_id=r.axis_id, matched_N=r.matched_N, null=st,
                             delta_obs=x, p=fn(x),
                             n_bg=cell_n[key], p_floor=cell_floor[key],
                             at_floor=bool(x > 0 and fn(x) <= cell_floor[key] + 1e-15)))
    sig = pd.DataFrame(recs)
    if sig.empty:
        sig["q"] = []
        sig["survives"] = []
        return sig
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


def run(args):
    metric, prefix, lane_name = LANES[args.lane]
    if args.faa:
        fc = fosmid_orf_counts_from_faa(args.faa)
    elif args.evidence:
        fc = fosmid_orf_counts_from_evidence(args.evidence)
    else:
        raise SystemExit("need --faa or --evidence for ORF counts")

    nul = load_nulls(args.nulls_dir, metric, prefix)
    obs = _load_obs(args.report, metric)
    obs["n_orfs"] = obs["fosmid"].astype(str).map(fc)
    obs = obs[obs["n_orfs"].notna()].copy()
    obs["matched_N"] = obs["n_orfs"].astype(int).map(nearest_n)
    print(f"[sig:{prefix}] {obs['fosmid'].nunique()} fosmids x "
          f"{obs['axis_id'].nunique()} axes ({len(obs)} obs rows)", flush=True)

    sig = score(nul, obs, metric)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(args.out, sep="\t", index=False)
    if not sig.empty:
        nfloor = int(sig[(sig.delta_obs > 0)]["at_floor"].sum())
        ntested = int((sig.delta_obs > 0).sum())
        nsurv = int(sig["survives"].sum())
        print(f"[sig:{prefix}] survivor cells q<{Q_THRESHOLD}: {nsurv}; "
              f"{nfloor}/{ntested} positive-delta at floor (ranked by effect size)",
              flush=True)
    print(f"[sig:{prefix}] wrote {args.out}", flush=True)


# --- selftest: reproduce scadc sig_emp_{reff,ieff}.tsv exactly --------------

def selftest(args):
    scadc = Path(args.scadc_root)
    cache = scadc / "cache"
    ref_dir = scadc / "reff"
    faa = Path("/home/tony/agentic_workspace/data/scadc/gene_centric/"
               "amino_acid/fosmids.spades_meta.gt29kb.faa")
    fc = fosmid_orf_counts_from_faa(faa)
    worst_overall = 0.0
    for lane in ("reff", "ieff"):
        metric, prefix, _ = LANES[lane]
        nul = load_nulls(cache, metric, prefix)
        obs = _load_obs(cache / f"{prefix}_axes_report.tsv", metric)
        obs["n_orfs"] = obs["fosmid"].astype(str).map(fc)
        obs = obs[obs["n_orfs"].notna()].copy()
        obs["matched_N"] = obs["n_orfs"].astype(int).map(nearest_n)
        sig = score(nul, obs, metric)
        ref = pd.read_csv(ref_dir / f"sig_emp_{prefix}.tsv", sep="\t")
        # join on the cell key, compare p / q / survives
        keys = ["fosmid", "element", "axis_id", "null"]
        m = sig.merge(ref, on=keys, suffixes=("_mine", "_ref"))
        dp = float(np.abs(m["p_mine"] - m["p_ref"]).max()) if len(m) else np.nan
        dq = float(np.abs(m["q_mine"] - m["q_ref"]).max()) if len(m) else np.nan
        surv_agree = int((m["survives_mine"] == m["survives_ref"]).sum())
        worst_overall = max(worst_overall, dp if np.isfinite(dp) else 0.0,
                            dq if np.isfinite(dq) else 0.0)
        print(f"[selftest:{prefix}] rows mine={len(sig)} ref={len(ref)} "
              f"matched={len(m)}; max|dp|={dp:.2e} max|dq|={dq:.2e}; "
              f"survives agree {surv_agree}/{len(m)}", flush=True)
    print(f"[selftest] worst |d(p or q)| overall = {worst_overall:.2e}")
    if worst_overall > 1e-9:
        raise SystemExit(f"selftest parity FAILED: {worst_overall:.2e}")
    print("[selftest] PASS")


def parse_args():
    ap = argparse.ArgumentParser(description="ECSPr significance (empirical-CCDF vs reused null)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="score an observed report vs the frozen null")
    r.add_argument("--lane", choices=list(LANES), required=True)
    r.add_argument("--report", required=True, help="ecspr_network solve output tsv")
    r.add_argument("--nulls-dir", required=True, help="dir with {reff,ieff}_null_canonical_N*.tsv")
    r.add_argument("--faa", default=None, help="prodigal proteins (>contig_N) for ORF counts")
    r.add_argument("--evidence", default=None, help="evidence parquet fallback for ORF counts")
    r.add_argument("--out", required=True)
    r.set_defaults(func=run)

    s = sub.add_parser("selftest", help="reproduce scadc sig_emp_{reff,ieff}.tsv")
    s.add_argument("--scadc-root", default="/home/tony/agentic_workspace/projects/scadc/"
                   "metabolic-modelling/main/metabolic-modelling/04_reaction_network")
    s.set_defaults(func=selftest)
    return ap.parse_args()


if __name__ == "__main__":
    a = parse_args()
    a.func(a)
