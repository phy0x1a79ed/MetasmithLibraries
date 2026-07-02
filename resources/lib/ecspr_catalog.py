"""ECSPr five-parquet reaction catalog -- MetaNetX x evidence synthesis.

Self-contained builder for the ECSPr reaction catalog the fabfos metasmith
pipeline consumes. Unlike `ecspr_solver.py` (a port of existing scadc solver
modules), this catalog has NO scadc source -- it is SYNTHESIZED here from three
inputs:

  * MetaNetX MNXref `reac_prop.tsv` / `chem_prop.tsv` (reaction equations,
    EC classifs, balance/transport flags; metabolite name/formula/mass/inchikey).
  * The stage-1 `evidence_table` parquet (per-ORF x per-channel MNXR nominations;
    03_layer2_evidence output).
  * The `evidence_weights` parquet (per (source, mnxr) belief weights E_full /
    E_dlec / n_orf; 04_reaction_network/10_evidence_weights.py output). That
    script splits each ORF's unit belief equally across its annotating lanes, then
    by raw_score within a lane and evenly across each nomination's MNXR fanout
    (conservation: per-ORF contributions sum to 1.0). `reaction_redundancy` is the
    per-reaction view of that allocation, augmented with the raw multiplicity
    (distinct channels / distinct orf x channel) the evidence table carries.

Five parquets are emitted into <out>, all join-able against an ORF set to subset:
  1. reactions.parquet           mnxr, ec_list, is_balanced, is_transport, equation
  2. metabolites.parquet         mnxm, name, formula, charge, mass, inchikey
  3. rxn_edges.parquet           mnxr, mnxm, coeff, side, compartment  (bipartite)
  4. orf_reactions.parquet       source, orf, mnxr, channel, intermediate_id,
                                 raw_score, projection_via  (NOT deduped by mnxr)
  5. reaction_redundancy.parquet source, mnxr, E_full, E_dlec, n_orf,
                                 n_channels, n_orf_channel  (one row / source,mnxr)

Equation grammar (reac_prop `mnx_equation`), parsed by parse_equation():
    "1 MNXM01@MNXD1 + 2 WATER@MNXD1 = 1 MNXM3@MNXD1"
    - ` = ` splits substrates (left) from products (right)
    - ` + ` splits terms on a side
    - a term is `<coeff> <metabolite>@<compartment>`; coeff may be a float; the
      metabolite id is the token bearing `@` (MNXM.../WATER/BIOMASS/...); the
      compartment is the `@`-suffix (MNXD1, MNXC.., ...). EMPTY / blank sides skip.

Env: pandas + pyarrow + stdlib only (no torch/networkx). Runs under e.g.
`scadc-metabolic-model`. chem_prop.tsv is ~810MB -- it is read in filtered
chunks against the metabolite set the parsed edges actually need.
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import pandas as pd

# reac_prop / chem_prop share a `#`-comment banner; the `#ID...` header line is
# itself a comment, so we skip all `#` lines and supply names explicitly.
REAC_PROP_COLS = ["mnxr", "mnx_equation", "reference", "classifs",
                  "is_balanced", "is_transport"]
CHEM_PROP_COLS = ["mnxm", "name", "reference", "formula", "charge", "mass",
                  "InChI", "InChIKey", "SMILES"]
# MetaNetX balance code: 'B' == mass/charge balanced. Transport code: 'T'.
BALANCED_CODE = "B"
TRANSPORT_CODE = "T"
# Filtered read of the 810MB chem_prop.tsv.
CHEM_CHUNK = 200_000


# =====================================================================
# Equation parsing
# =====================================================================

def _parse_term(term: str):
    """'2 WATER@MNXD1' -> (mnxm, coeff, compartment) or None for blanks.

    The metabolite token is the one bearing '@'; everything before it is the
    coefficient (a single float in practice, but joined defensively to tolerate a
    multi-token coefficient). Returns coeff as float, NaN if unparseable."""
    toks = term.split()
    if not toks:
        return None
    met_tok = None
    met_pos = None
    for i, t in enumerate(toks):
        if "@" in t:
            met_tok = t
            met_pos = i
            break
    if met_tok is None:
        return None
    mnxm, _, compartment = met_tok.partition("@")
    if not mnxm or mnxm == "EMPTY":
        return None
    coeff_txt = " ".join(toks[:met_pos]).strip()
    try:
        coeff = float(coeff_txt) if coeff_txt else 1.0
    except ValueError:
        coeff = float("nan")
    return mnxm, coeff, compartment


def parse_equation(mnxr: str, equation):
    """Yield (mnxr, mnxm, coeff, side, compartment) rows for one equation.

    side in {'substrate','product'} (left/right of ' = '). Handles EMPTY / NaN /
    blank sides and the '@MNXDn' compartment suffix."""
    if not isinstance(equation, str) or "=" not in equation:
        return
    left, _, right = equation.partition("=")
    for side, chunk in (("substrate", left), ("product", right)):
        for term in chunk.split(" + "):
            parsed = _parse_term(term)
            if parsed is None:
                continue
            mnxm, coeff, compartment = parsed
            yield (mnxr, mnxm, coeff, side, compartment)


# =====================================================================
# Loaders
# =====================================================================

def load_reac_prop(path: Path, wanted_mnxr: set | None) -> pd.DataFrame:
    """reac_prop.tsv filtered to wanted_mnxr (all reactions if None)."""
    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     names=REAC_PROP_COLS, dtype=str, na_filter=False)
    if wanted_mnxr is not None:
        df = df[df["mnxr"].isin(wanted_mnxr)]
    return df.reset_index(drop=True)


def load_chem_prop(path: Path, wanted_mnxm: set) -> pd.DataFrame:
    """chem_prop.tsv read in chunks, kept only where mnxm in wanted_mnxm.

    usecols drops the heavy InChI/SMILES columns up front; the row filter keeps
    memory bounded on the 810MB file."""
    keep = ["mnxm", "name", "formula", "charge", "mass", "InChIKey"]
    parts = []
    reader = pd.read_csv(path, sep="\t", comment="#", header=None,
                         names=CHEM_PROP_COLS, usecols=keep, dtype=str,
                         na_filter=False, chunksize=CHEM_CHUNK)
    for chunk in reader:
        hit = chunk[chunk["mnxm"].isin(wanted_mnxm)]
        if len(hit):
            parts.append(hit)
    if parts:
        out = pd.concat(parts, ignore_index=True)
    else:
        out = pd.DataFrame(columns=keep)
    out = out.rename(columns={"InChIKey": "inchikey"})
    out["charge"] = pd.to_numeric(out["charge"], errors="coerce")
    out["mass"] = pd.to_numeric(out["mass"], errors="coerce")
    return out[["mnxm", "name", "formula", "charge", "mass", "inchikey"]]


# =====================================================================
# Build
# =====================================================================

def build(evidence_path: Path, weights_path: Path, reac_prop_path: Path,
          chem_prop_path: Path, out_dir: Path,
          reaction_set: set | None = None) -> dict:
    """Synthesize the five parquets into out_dir. Returns {name: row_count}."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ev = pd.read_parquet(evidence_path)
    weights = pd.read_parquet(weights_path)

    # reaction universe: evidence MNXR unless an explicit set is passed.
    wanted_mnxr = set(reaction_set) if reaction_set is not None \
        else set(ev["mnxr"].unique())

    # ---- reac_prop -> equations + rxn_edges -----------------------------
    reac = load_reac_prop(reac_prop_path, wanted_mnxr)
    edge_rows = []
    for mnxr, eq in zip(reac["mnxr"], reac["mnx_equation"]):
        edge_rows.extend(parse_equation(mnxr, eq))
    rxn_edges = pd.DataFrame(
        edge_rows, columns=["mnxr", "mnxm", "coeff", "side", "compartment"])

    # (1) reactions: keep only reactions that produced >=1 parsed edge, so every
    # reactions row is guaranteed present in rxn_edges (the catalog invariant).
    with_edges = set(rxn_edges["mnxr"].unique())
    reactions = reac[reac["mnxr"].isin(with_edges)].copy()
    reactions["ec_list"] = reactions["classifs"]
    reactions["is_balanced"] = reactions["is_balanced"].eq(BALANCED_CODE)
    reactions["is_transport"] = reactions["is_transport"].eq(TRANSPORT_CODE)
    reactions = reactions.rename(columns={"mnx_equation": "equation"})
    reactions = reactions[["mnxr", "ec_list", "is_balanced",
                           "is_transport", "equation"]].reset_index(drop=True)
    rxn_edges = rxn_edges[rxn_edges["mnxr"].isin(with_edges)].reset_index(drop=True)

    # (2) metabolites: only MNXM that appear in the kept edges.
    wanted_mnxm = set(rxn_edges["mnxm"].unique())
    metabolites = load_chem_prop(chem_prop_path, wanted_mnxm)

    # (3) rxn_edges already built above.

    # (4) orf_reactions: raw multiplicity, NOT deduped by mnxr.
    orf_cols = ["source", "orf", "mnxr", "channel", "intermediate_id",
                "raw_score", "projection_via"]
    orf_reactions = ev[orf_cols].copy()
    if reaction_set is not None:
        orf_reactions = orf_reactions[orf_reactions["mnxr"].isin(wanted_mnxr)]
    orf_reactions = orf_reactions.reset_index(drop=True)

    # (5) reaction_redundancy: evidence-weights view + raw redundancy counts.
    ev_scope = ev[ev["mnxr"].isin(wanted_mnxr)] if reaction_set is not None else ev
    n_channels = (ev_scope.groupby(["source", "mnxr"])["channel"]
                  .nunique().rename("n_channels").reset_index())
    oc = ev_scope.drop_duplicates(["source", "mnxr", "orf", "channel"])
    n_orf_channel = (oc.groupby(["source", "mnxr"]).size()
                     .rename("n_orf_channel").reset_index())
    w = weights
    if reaction_set is not None:
        w = w[w["mnxr"].isin(wanted_mnxr)]
    redundancy = (w.merge(n_channels, on=["source", "mnxr"], how="outer")
                  .merge(n_orf_channel, on=["source", "mnxr"], how="outer"))
    for c in ("E_full", "E_dlec"):
        redundancy[c] = redundancy[c].fillna(0.0)
    for c in ("n_orf", "n_channels", "n_orf_channel"):
        redundancy[c] = redundancy[c].fillna(0).astype(int)
    redundancy = redundancy[["source", "mnxr", "E_full", "E_dlec", "n_orf",
                             "n_channels", "n_orf_channel"]].reset_index(drop=True)

    outputs = {
        "reactions": reactions,
        "metabolites": metabolites,
        "rxn_edges": rxn_edges,
        "orf_reactions": orf_reactions,
        "reaction_redundancy": redundancy,
    }
    counts = {}
    for name, df in outputs.items():
        df.to_parquet(out_dir / f"{name}.parquet", index=False)
        counts[name] = len(df)
    return counts


# =====================================================================
# Canonical paths (self-test fixtures)
# =====================================================================

_SCADC = Path("/home/tony/agentic_workspace")
CANON_EVIDENCE = _SCADC / ("projects/scadc/metabolic-modelling/main/"
                           "metabolic-modelling/03_layer2_evidence/cache/"
                           "evidence_table_dlec.parquet")
CANON_WEIGHTS = _SCADC / ("projects/scadc/metabolic-modelling/main/"
                          "metabolic-modelling/04_reaction_network/cache/"
                          "evidence_weights.parquet")
CANON_REAC_PROP = _SCADC / "data/scadc/references/metanetx/reac_prop.tsv"
CANON_CHEM_PROP = _SCADC / "data/scadc/references/metanetx/chem_prop.tsv"


def _selftest() -> int:
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        counts = build(CANON_EVIDENCE, CANON_WEIGHTS, CANON_REAC_PROP,
                       CANON_CHEM_PROP, out)
        print("row counts:")
        for name, n in counts.items():
            print(f"  {name:22s} {n:,}")

        reactions = pd.read_parquet(out / "reactions.parquet")
        metabolites = pd.read_parquet(out / "metabolites.parquet")
        rxn_edges = pd.read_parquet(out / "rxn_edges.parquet")
        orf_reactions = pd.read_parquet(out / "orf_reactions.parquet")
        redundancy = pd.read_parquet(out / "reaction_redundancy.parquet")

        for name, df in (("reactions", reactions), ("metabolites", metabolites),
                         ("rxn_edges", rxn_edges),
                         ("orf_reactions", orf_reactions),
                         ("reaction_redundancy", redundancy)):
            assert len(df) > 0, f"{name} is empty"

        # every reactions mnxr appears in rxn_edges
        r_set = set(reactions["mnxr"])
        e_set = set(rxn_edges["mnxr"])
        missing = r_set - e_set
        assert not missing, f"{len(missing)} reactions absent from rxn_edges"

        # multiplicity preserved: orf_reactions == evidence table row count
        ev_n = len(pd.read_parquet(CANON_EVIDENCE))
        assert len(orf_reactions) == ev_n, \
            f"orf_reactions {len(orf_reactions)} != evidence {ev_n}"

        # one row per (source, mnxr) in reaction_redundancy
        dup = redundancy.duplicated(["source", "mnxr"]).sum()
        assert dup == 0, f"reaction_redundancy has {dup} dup (source,mnxr)"

        # metabolites: every edge mnxm should mostly resolve (report coverage)
        cov = len(set(metabolites["mnxm"]) & set(rxn_edges["mnxm"]))
        print(f"metabolite coverage: {cov:,} / "
              f"{rxn_edges['mnxm'].nunique():,} distinct edge metabolites")
    print("\nPASS")
    return 0


# =====================================================================
# CLI
# =====================================================================

def _cli(argv) -> int:
    ap = argparse.ArgumentParser(description="ECSPr five-parquet reaction catalog")
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="synthesize the five parquets")
    b.add_argument("--evidence", required=True)
    b.add_argument("--weights", required=True)
    b.add_argument("--reac-prop", required=True)
    b.add_argument("--chem-prop", required=True)
    b.add_argument("--out", required=True)
    # cross-ref tables are optional -- only needed if cross-ref columns are added.
    b.add_argument("--reac-xref", default=None)
    b.add_argument("--chem-xref", default=None)

    sub.add_parser("selftest", help="build on canonical fixtures + assert")

    args = ap.parse_args(argv)
    if args.cmd == "selftest":
        return _selftest()
    if args.cmd == "build":
        counts = build(Path(args.evidence), Path(args.weights),
                       Path(args.reac_prop), Path(args.chem_prop), Path(args.out))
        for name, n in counts.items():
            print(f"{name:22s} {n:,}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(_cli(sys.argv[1:]))
