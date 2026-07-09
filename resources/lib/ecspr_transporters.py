"""ECSPr transporter subsystem -- TC-family/number -> substrate -> MNXM resolution.

Turns a TCDB-homology transporter annotation into the per-organism
{organism: {MNXM: belief}} map the community graph consumes to gate membrane
crossing. This is the recall/precision knob of the whole method, kept behind a clean
API so a coarser or richer substrate mapping swaps in without touching the graph.

Substrate resolution (the crux):
  tc_number --(TCDB substrate table)--> ChEBI ids --(MetaNetX chem_xref)--> MNXM
  with a family-level fallback: if a hit's exact 5-field TC number has no substrate
  record, union the substrates of every TC number sharing its 3-field family. This is
  recall-oriented (the user's explicit goal) and flagged per row via `resolved_via`.

Cross-check lane (optional): the MetaNetX is_transport reactions the organism's
reaction annotation nominates directly name the metabolite that crosses; unioning
those MNXM in raises recall for substrates TCDB's substrate table misses. Wired via
`metanetx_transport_substrates`.

Belief: each transporter ORF hit contributes evidence e in (0,1] (default pident/100);
per (organism, MNXM) the belief is the noisy-OR 1 - prod(1 - e) over supporting hits
-- bounded in [0,1), monotonic in evidence, saturating with transporter copy number
(more transporter genes -> more transport capacity, but never runaway conductance).

Env: pandas + numpy (CPU).
"""
from __future__ import annotations

import argparse
import pickle
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

_CHEBI_RE = re.compile(r"chebi:(\d+)", re.IGNORECASE)


def tc_family(tc_number: str) -> str:
    """5-field TC number -> 3-field family, e.g. 3.A.1.104.2 -> 3.A.1."""
    return ".".join(str(tc_number).split(".")[:3])


# =====================================================================
# reference crosswalks
# =====================================================================

def load_chebi_to_mnxm(chem_xref: Path) -> dict:
    """{chebi_number(str): set(MNXM)} from MetaNetX chem_xref.tsv (source\tMNXM\tdesc).
    ChEBI appears as CHEBI:/chebi: -- matched case-insensitively on the number."""
    out: dict = defaultdict(set)
    with open(chem_xref) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, mnxm = parts[0], parts[1]
            if not mnxm.startswith("MNXM"):
                continue
            m = _CHEBI_RE.match(src)
            if m:
                out[m.group(1)].add(mnxm)
    return {k: v for k, v in out.items()}


def load_tc_substrates(tcdb_substrates: Path) -> dict:
    """{tc_number: set(chebi_number)} from the TCDB substrate dump
    (tc_number \\t CHEBI:id;label|CHEBI:id;label...)."""
    out: dict = defaultdict(set)
    with open(tcdb_substrates) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            tc, subs = parts[0], parts[1]
            for entry in subs.split("|"):
                m = _CHEBI_RE.search(entry)
                if m:
                    out[tc].add(m.group(1))
    return {k: v for k, v in out.items()}


def build_tc_to_mnxm(tc_substrates: dict, chebi2mnxm: dict) -> tuple[dict, dict]:
    """Return (tc_number -> set(MNXM), family -> set(MNXM)). The family map unions
    every member TC number's substrates for the recall-oriented fallback."""
    tc2mnxm: dict = {}
    fam2mnxm: dict = defaultdict(set)
    for tc, chebis in tc_substrates.items():
        mnxms = set()
        for c in chebis:
            mnxms |= chebi2mnxm.get(c, set())
        if mnxms:
            tc2mnxm[tc] = mnxms
            fam2mnxm[tc_family(tc)] |= mnxms
    return tc2mnxm, dict(fam2mnxm)


# =====================================================================
# per-organism transporter belief map
# =====================================================================

def resolve_hit_mnxm(tc_number: str, tc2mnxm: dict, fam2mnxm: dict) -> tuple[set, str]:
    """(MNXM set, resolved_via) for one TC number: exact number first, then family."""
    if tc_number in tc2mnxm:
        return tc2mnxm[tc_number], "tc_number"
    fam = tc_family(tc_number)
    if fam in fam2mnxm:
        return fam2mnxm[fam], "tc_family"
    return set(), "unresolved"


def build_transporter_map(hits: pd.DataFrame, tc2mnxm: dict, fam2mnxm: dict,
                          score_col: str = "pident", score_scale: float = 100.0,
                          min_evidence: float = 0.25,
                          metanetx_transport: dict | None = None) -> tuple[dict, pd.DataFrame]:
    """{organism: {MNXM: belief}} via noisy-OR over transporter ORF hits.

    hits: tcdb_hits parquet (organism, orf, tc_number, pident, ...).
    metanetx_transport: optional {organism: set(MNXM)} cross-check lane (belief floor
        `min_evidence` for those, unioned in).
    Returns (belief_map, audit_df).
    """
    # accumulate per (org, mnxm): list of evidence e in (0,1]
    ev: dict = defaultdict(lambda: defaultdict(list))
    audit_rows = []
    for r in hits.itertuples(index=False):
        e = float(getattr(r, score_col)) / score_scale
        e = min(max(e, 0.0), 0.999)
        mnxms, via = resolve_hit_mnxm(r.tc_number, tc2mnxm, fam2mnxm)
        audit_rows.append({"organism": r.organism, "orf": r.orf,
                           "tc_number": r.tc_number, "n_mnxm": len(mnxms),
                           "resolved_via": via, "evidence": e})
        for m in mnxms:
            ev[r.organism][m].append(e)

    belief: dict = defaultdict(dict)
    for org, mm in ev.items():
        for m, es in mm.items():
            noisy_or = 1.0 - np.prod([1.0 - e for e in es])
            belief[org][m] = float(noisy_or)

    # optional MetaNetX is_transport cross-check lane
    if metanetx_transport:
        for org, mnxms in metanetx_transport.items():
            for m in mnxms:
                cur = belief[org].get(m, 0.0)
                belief[org][m] = float(1.0 - (1.0 - cur) * (1.0 - min_evidence))

    return dict(belief), pd.DataFrame(audit_rows)


# =====================================================================
# MetaNetX is_transport substrate lane (optional cross-check)
# =====================================================================

def load_transport_reaction_substrates(reac_prop: Path) -> dict:
    """{MNXR: set(MNXM)} for is_transport reactions -- the metabolite(s) that cross.
    Parses the compartment-tagged equation; the crossing species is any MNXM that
    appears on both sides in different @MNXD compartments (i.e. is transported)."""
    out: dict = {}
    eq_re = re.compile(r"(MNXM\w+)@(MNXD\d+)")
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[5] != "T":
                continue
            mnxr, eq = parts[0], parts[1]
            comp_of = defaultdict(set)
            for mnxm, comp in eq_re.findall(eq):
                comp_of[mnxm].add(comp)
            crossing = {m for m, cs in comp_of.items() if len(cs) > 1}
            if crossing:
                out[mnxr] = crossing
    return out


# =====================================================================
# CLI
# =====================================================================

def cmd_build(args):
    print("[transporters] loading crosswalks...", flush=True)
    chebi2mnxm = load_chebi_to_mnxm(Path(args.chem_xref))
    tc_subs = load_tc_substrates(Path(args.tc_substrates))
    tc2mnxm, fam2mnxm = build_tc_to_mnxm(tc_subs, chebi2mnxm)
    print(f"  ChEBI->MNXM: {len(chebi2mnxm):,} ChEBI  |  "
          f"TC substrate records: {len(tc_subs):,}  |  "
          f"TC#->MNXM: {len(tc2mnxm):,}  families: {len(fam2mnxm):,}", flush=True)

    hits = pd.read_parquet(args.hits)
    mnx_transport = None
    if args.reac_prop and args.org_reactions:
        rxn_subs = load_transport_reaction_substrates(Path(args.reac_prop))
        org_rxn = pickle.load(open(args.org_reactions, "rb"))  # {org: {mnxr: E}}
        mnx_transport = {}
        for org, w in org_rxn.items():
            s = set()
            for mnxr in w:
                s |= rxn_subs.get(mnxr, set())
            mnx_transport[org] = s
        print("  MetaNetX is_transport lane: "
              + ", ".join(f"{o}={len(s)}" for o, s in mnx_transport.items()), flush=True)

    belief, audit = build_transporter_map(hits, tc2mnxm, fam2mnxm,
                                          metanetx_transport=mnx_transport)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wb") as fh:
        pickle.dump(belief, fh, protocol=pickle.HIGHEST_PROTOCOL)
    if args.audit:
        audit.to_parquet(args.audit, index=False)

    print("\n[transporters] per-organism substrate coverage:")
    for org in sorted(belief):
        print(f"  {org}: {len(belief[org]):,} transportable MNXM  "
              f"(mean belief {np.mean(list(belief[org].values())):.3f})")
    via = audit["resolved_via"].value_counts().to_dict()
    print(f"  resolution: {via}")
    print(f"[transporters] wrote {args.out}")


def parse_args():
    ap = argparse.ArgumentParser(description="TC -> substrate -> MNXM transporter map")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build"); p.set_defaults(fn=cmd_build)
    p.add_argument("--hits", required=True, help="tcdb_hits.parquet")
    p.add_argument("--tc-substrates", required=True, help="tcdb_substrates.tsv")
    p.add_argument("--chem-xref", required=True, help="MetaNetX chem_xref.tsv")
    p.add_argument("--reac-prop", default=None, help="MetaNetX reac_prop.tsv (is_transport lane)")
    p.add_argument("--org-reactions", default=None, help="pickle {org:{mnxr:E}} for the is_transport lane")
    p.add_argument("--out", required=True, help="output pickle {org:{mnxm:belief}}")
    p.add_argument("--audit", default=None, help="optional audit parquet")
    return ap.parse_args()


if __name__ == "__main__":
    a = parse_args()
    a.fn(a)
