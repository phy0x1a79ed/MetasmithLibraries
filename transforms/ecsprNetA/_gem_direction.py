"""Network A per-reaction direction ratios from the curated GEM's NATIVE flux bounds.

Network B and the frozen reference derive per-reaction directionality from a
THERMODYNAMIC ensemble: a continuous ratio ``g_rev/g_fwd = exp(dG'/RT)`` fused from
eQuilibrator + dGbyG + BioCyc. Network A is built from a curated genome-scale model
(iECDH10B) whose reactions ALREADY encode direction natively, as flux bounds -- so the
honest, native directionality for Network A is those bounds, no thermodynamics and no
external ensemble.

BOUNDS -> RATIO
---------------
``ratio`` is ``g_rev/g_fwd`` where "fwd" is the substrate->product orientation the
solver's MetaNetX roles define (see ``ecspr_directed.orient_and_weight``: substrate
edge met->rxn is the tail, product edge rxn->met is the head, ``gm = ratio * gp``).
A reaction's permitted-direction set is read from its bounds -- forward feasible iff
``ub > 0``, reverse feasible iff ``lb < 0`` -- and mapped:

    reversible           lb < 0 < ub      ratio = 1.0        symmetric: gm == gp, the
                                                             proven no-op == undirected
    forward-irreversible ub > 0, lb >= 0  ratio = FLOOR      strong forward diode
    reverse-irreversible ub <= 0, lb < 0  ratio = 1/FLOOR    strong reverse diode
    blocked              lb == ub == 0     (row omitted)     the solver defaults a
                                                             missing MNXR to ratio 1.0

FLOOR is the engine's ``DIODE_BACKWARD_FLOOR`` (imported, never restated): the sharpest
backward throttle the rectified-network primitive honours before its OWN safety clip
(``directed_ceff`` clips ``gm`` up to ``FLOOR * gp`` to keep the grounded Laplacian
nonsingular). A forward-irreversible ratio of exactly FLOOR therefore lands on that
clip -- the strongest one-way diode the solver represents without a singular system --
and its reverse-irreversible reciprocal ``1/FLOOR`` is the mirror on the other branch.

This is a BINARY model (reversible / irreversible), unlike Network B's CONTINUOUS
``exp(dG'/RT)``: a GEM bound states only the SIGN of feasible flux, not a free-energy
magnitude, so the diode has ONE fixed magnitude rather than a per-reaction one. The
consequence is extreme ratios (``FLOOR`` / ``1/FLOOR`` ~ 1e-9 / 1e9) -- exactly the
backflow regime the directed solver is being hardened for.

Duplicate GEM reactions crosswalking to one MNXR are combined by the UNION of their
permitted directions (a reaction is forward-feasible for the MNXR if ANY copy is, and
likewise reverse). So a reversible copy beside an irreversible one, or two copies
irreversible in opposite senses, jointly permit both directions -> reversible (1.0):
direction is only imposed where every GEM copy agrees it is one-way.

ORIENTATION CAVEAT
------------------
The GEM's ``lb``/``ub`` are in the GEM reaction's OWN equation orientation. This builder
assumes that orientation matches the MetaNetX canonical (``reac_prop``) orientation the
solver's roles use. Network B's annotator explicitly ALIGNS BioCyc direction to MNXR;
this builder does not re-align GEM bounds. Where a crosswalked reaction's BiGG equation
is written reverse to its MNXR equation, its irreversible ratio would be inverted. This
is the one modelling simplification of the native-bounds lane vs the aligned ensemble;
it is flagged for the lead's reconciliation, not silently assumed away.

Canon-free: every path is a CLI argument; no SCADC path and no ``canon`` import. The
rxn->MNXR crosswalk is REUSED verbatim from ``ecspr_network.gem_crosswalk`` so the
emitted MNXR keys are byte-identical to the Network A base graph's ``("rxn", mnxr)``
nodes -- a direction row keyed on an MNXR the base never carries would silently be a
no-op. The output is exactly the two columns ``ecspr_network.load_direction_ratios``
reads: ``mnxr`` (str), ``ratio`` (float).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# resources/lib is <repo>/resources/lib; this file is <repo>/transforms/ecsprNetA/.
# Add it so the crosswalk and the diode floor are REUSED, not re-derived: identical
# MNXR keys as the base builder, and the diode magnitude as one named engine constant.
_LIB = Path(__file__).resolve().parents[2] / "resources" / "lib"
sys.path.insert(0, str(_LIB))
from ecspr_network import gem_crosswalk            # noqa: E402  rxn_id -> canonical MNXR
from ecspr_directed import DIODE_BACKWARD_FLOOR    # noqa: E402  the diode magnitude

ELEMENTS = ["C", "N", "S", "P"]


def _rxn_bounds(model_json: Path) -> dict:
    """{rxn_id: (lower_bound, upper_bound)} from the curated GEM JSON."""
    import cobra                                   # lazy: only this path needs it
    m = cobra.io.load_json_model(str(model_json))
    return {r.id: (float(r.lower_bound), float(r.upper_bound)) for r in m.reactions}


def _feasible(lb: float, ub: float) -> tuple[bool, bool]:
    """(forward_feasible, reverse_feasible) from flux bounds. Forward flux is possible
    iff the upper bound is positive; reverse iff the lower bound is negative."""
    return (ub > 0.0, lb < 0.0)


def build_direction(model_json, reac_xref, bipartite_dir, elements=ELEMENTS,
                    floor: float = DIODE_BACKWARD_FLOOR):
    """Emit the (mnxr, ratio) direction table for Network A from GEM bounds.

    Returns (DataFrame[mnxr, ratio], stats). The crosswalk (and its universe-preference)
    is identical to the base builder's, so every emitted MNXR is a real Network A node.
    """
    xw, stats = gem_crosswalk(Path(model_json), Path(reac_xref), Path(bipartite_dir),
                              elements)
    bounds = _rxn_bounds(Path(model_json))

    # Per-MNXR union of permitted directions across every GEM reaction mapping to it.
    fwd: dict[str, bool] = {}
    rev: dict[str, bool] = {}
    for row in xw.itertuples():
        lb, ub = bounds.get(row.rxn_id, (0.0, 0.0))
        f, r = _feasible(lb, ub)
        fwd[row.mnxr] = fwd.get(row.mnxr, False) or f
        rev[row.mnxr] = rev.get(row.mnxr, False) or r

    rows = []
    for mnxr in sorted(fwd):
        f, r = fwd[mnxr], rev[mnxr]
        if f and r:
            ratio = 1.0                    # reversible -> symmetric no-op
        elif f:
            ratio = float(floor)           # forward-irreversible -> strong forward diode
        elif r:
            ratio = 1.0 / float(floor)     # reverse-irreversible -> strong reverse diode
        else:
            continue                       # blocked -> omit -> solver default reversible
        rows.append((mnxr, ratio))

    df = pd.DataFrame(rows, columns=["mnxr", "ratio"])
    return df, stats


def _validate(df: pd.DataFrame) -> None:
    """The load_direction_ratios contract: exactly (mnxr, ratio), unique keys, ratio a
    strictly-positive conductance ratio. A self-check so a malformed table fails here,
    at build time, not deep inside the directed solve."""
    assert list(df.columns) == ["mnxr", "ratio"], f"columns are {list(df.columns)}"
    assert df["mnxr"].is_unique, "duplicate MNXR keys -- the union collapse failed"
    assert (df["ratio"] > 0).all(), "ratio must be strictly positive (a conductance ratio)"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Network A per-reaction direction ratios from curated-GEM flux bounds")
    ap.add_argument("--model", required=True, help="curated GEM JSON (iECDH10B)")
    ap.add_argument("--reac-xref", required=True, help="MetaNetX reac_xref.tsv (current)")
    ap.add_argument("--bipartite-dir", required=True,
                    help="per-element mnx_bipartite_{X}.pkl dir (crosswalk universe pref)")
    ap.add_argument("--out", required=True, help="output parquet (columns: mnxr, ratio)")
    ap.add_argument("--elements", nargs="+", default=ELEMENTS)
    ap.add_argument("--floor", type=float, default=DIODE_BACKWARD_FLOOR,
                    help="diode magnitude; default the engine DIODE_BACKWARD_FLOOR")
    a = ap.parse_args()

    df, stats = build_direction(a.model, a.reac_xref, a.bipartite_dir, a.elements, a.floor)
    _validate(df)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)

    n_fwd = int((df["ratio"] < 1.0).sum())
    n_rev = int((df["ratio"] > 1.0).sum())
    n_sym = int((df["ratio"] == 1.0).sum())
    print(f"[gem-direction] {stats['n_resolved']}/{stats['n_reactions']} GEM reactions "
          f"crosswalked -> {len(df)} unique MNXR -> {out}")
    print(f"[gem-direction] forward-irrev {n_fwd} (ratio={a.floor:g}); "
          f"reverse-irrev {n_rev} (ratio={1.0 / a.floor:g}); "
          f"reversible {n_sym} (ratio=1.0); blocked omitted (solver default reversible)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
