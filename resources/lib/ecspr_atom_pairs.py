"""Atom-transfer PAIRS from cached atom-mapping output: which metabolite's atoms
land in which metabolite, per element.

WHY THIS EXISTS
---------------
The bipartite met-rxn graph the solve runs on is a STAR: a reaction is a node joined
to every metabolite it touches. Eliminating that node -- which is exactly what the
Woodbury update's `_build_delta` does -- leaves a CLIQUE over every participant, so
every pair of participants gets conductance whether or not one atom passes between
them. The edge weights are atom counts, so the method looks stoichiometry-aware; the
TOPOLOGY throws the atom mapping away.

Measured on the real carbon graph: 52.8% of the reaction-induced metabolite pairs
transfer no atoms at all. On MNXR106432 (acetyl-CoA + CO2 + NADPH = pyruvate + NADP+
+ CoA) the pyruvate-NADPH channel carries ZERO carbons at conductance 0.700, while
the real one-carbon pyruvate-CO2 channel gets 0.033 -- a zero-atom channel 21x more
conductive than a real one.

The information needed to fix it already exists and is already computed. The unifier
(`20b_aam_unify.py:count_transits`) holds the substrate metabolite and the product
metabolite of each mapped atom IN THE SAME LOOP ITERATION, and then discards the
pairing into two independent role counters:

    sm = sub_mnxms[si]                       # substrate MNXM   <- the pair
    pm = prod_mnxms[j]                       # product MNXM     <- is right here
    if sm: counts[el][f"{sm}:substrate"] += 1     # and then thrown
    if pm: counts[el][f"{pm}:product"]   += 1     # away, separately

This module keeps `(sm, pm)`. Everything else is the same computation.

RE-EXTRACTION, NOT RE-MAPPING. The mapped SMILES are cached with their atom-map
numbers; nothing here runs RXNMapper or LocalMapper. This is arithmetic over a cache.

WHY A METABOLITE-LEVEL PAIR IS STILL NOT ENOUGH
-----------------------------------------------
Worth stating because it is the trap this whole design walks around. Pairing at the
METABOLITE level does not fix the stated problem. The true atom map splits
acetyl-CoA's 23 carbons by destination: atoms [0,1] -- the acetyl group -- become
pyruvate; the other 21 -- the CoA moiety -- become CoA. So acetyl-CoA -> CoA is a
GENUINE 21-carbon edge that any metabolite-level graph keeps, and a path arriving on
acetyl carbons could still leave through the CoA moiety carrying 21 carbons it never
brought. Only atom-resolved nodes deny that, which is why this module emits the atom
INDICES and not merely the counts -- the indices are what the class refinement
downstream partitions on.

WHAT IS REFUSED, AND WHY
------------------------
  * ATOM-UNBALANCED reactions. Generic metabolites with no SMILES were dropped before
    mapping, so the mapper re-routed their atoms onto whatever concrete product
    remained. Those reactions' COUNTS were presence-trustworthy -- which is why the
    incumbent could use them -- but their PAIRS are fabricated. Counts survive a lie
    about destination; pairs are nothing but destination.
  * DUPLICATE-METABOLITE reactions. `match_mols_to_mnxms` is greedy by canonical
    SMILES (`candidates.pop(0)`), so when one metabolite appears twice the
    substrate<->product assignment is arbitrary. Harmless at metabolite level;
    fabricates a split at atom level.

Both are counted and reported rather than dropped quietly -- the coverage cost is a
finding about the method's reach, not an embarrassment to hide.

Env: rdkit + pandas. No SCADC paths: every input arrives as an argument.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ELEMENTS = ("C", "N", "S", "P")
EQ_TERM = re.compile(r"(\d+(?:\.\d+)?)\s+(MNXM\d+|MNXM\w+|WATER|BIOMASS)")

PAIR_COLS = ("mnxr", "element", "substrate", "product", "n_atoms",
             "sub_idx", "prod_idx")
STATUS_COLS = ("mnxr", "status", "n_sub", "n_prod", "n_mapped_sub", "n_mapped_prod")


# =====================================================================
# MetaNetX side (same derivation the unifier uses)
# =====================================================================

def canon_smiles(smi: str):
    """Canonical SMILES with atom-map numbers stripped."""
    if not smi:
        return None
    try:
        mol = Chem.MolFromSmiles(smi)
    except Exception:
        return None
    if mol is None:
        return None
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def load_mnxm_smiles(chem_prop: Path, want: set) -> dict:
    out = {}
    with open(chem_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in want:
                continue
            if p[8].strip():
                out[p[0]] = p[8].strip()
    return out


def load_equations(reac_prop: Path, want: set) -> dict:
    out = {}
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 1 and p[0] in want:
                out[p[0]] = p[1]
    return out


def parse_equation(eq: str):
    """(substrates, products), each with multiplicity, in equation order."""
    if "=" not in eq:
        return None
    lhs, rhs = eq.split("=", 1)

    def expand(s):
        out = []
        for coef, m in EQ_TERM.findall(s):
            for _ in range(max(1, int(round(float(coef))))):
                out.append(m)
        return out

    L, R = expand(lhs), expand(rhs)
    return (L, R) if (L and R) else None


def match_mols(mols: list, mnxms: list, canon: dict):
    """Greedy canonical-SMILES match, aligned with `mols`.

    Also returns whether the match was AMBIGUOUS -- i.e. some canonical SMILES
    appears more than once among the MNXMs, so `pop(0)` picked arbitrarily. At
    metabolite level that is harmless; at atom level it invents a split.
    """
    pool = defaultdict(list)
    for m in mnxms:
        cs = canon.get(m)
        if cs:
            pool[cs].append(m)
    ambiguous = any(len(v) > 1 for v in pool.values())
    out = []
    for mol in mols:
        if mol is None:
            out.append(None)
            continue
        m2 = Chem.Mol(mol)
        for a in m2.GetAtoms():
            a.SetAtomMapNum(0)
        try:
            cs = Chem.MolToSmiles(m2, canonical=True)
        except Exception:
            out.append(None)
            continue
        cand = pool.get(cs)
        out.append(cand.pop(0) if cand else None)
    return out, ambiguous


# =====================================================================
# The pairing itself -- the two lines the unifier throws away
# =====================================================================

def canonical_ranks(mol):
    """Atom -> canonical rank, invariant to how this reaction happened to write it.

    THIS IS LOAD-BEARING, AND THE OBVIOUS THING IS WRONG. `a.GetIdx()` is the index
    within THIS reaction's template molecule, and the same metabolite is not written
    the same way in every reaction. Measured over the mapper universe: 248 of 847
    carbon-bearing metabolites that appear in more than one reaction have
    INCONSISTENT idx -> canonical-rank mappings. Pyruvate (MNXM23) alone shows up
    with 3 distinct atom orderings across 63 reactions.

    So `(metabolite, GetIdx())` is not an atom -- it is a different atom depending on
    which reaction you read it from, and a graph keyed on it silently welds unrelated
    atoms together and splits identical ones apart. `CanonicalRankAtoms` is invariant
    to input ordering by construction, so `(metabolite, rank)` IS an atom.
    """
    m2 = Chem.Mol(mol)
    for a in m2.GetAtoms():
        a.SetAtomMapNum(0)          # map numbers are per-reaction; they'd poison the rank
    try:
        return list(Chem.CanonicalRankAtoms(m2, breakTies=True))
    except Exception:
        return None


def pairs_from_mapped(mapped_smi: str, sub_mnxms: list, prod_mnxms: list,
                      canon: dict):
    """(pairs, status). pairs: {(element, sm, pm) -> [(sub_rank, prod_rank)]}.

    The atom identifiers are CANONICAL RANKS within each metabolite's own molecule,
    not RDKit template indices -- see `canonical_ranks`. That is what makes
    `(metabolite, rank)` a well-defined node across reactions, which is the whole
    premise of the atom graph.
    """
    if not mapped_smi or pd.isna(mapped_smi):
        return {}, "no_mapping"
    try:
        rxn = AllChem.ReactionFromSmarts(mapped_smi, useSmiles=True)
    except Exception:
        return {}, "unparseable"
    if rxn is None:
        return {}, "unparseable"

    n_sub_t = rxn.GetNumReactantTemplates()
    n_prod_t = rxn.GetNumProductTemplates()
    sub_mols = [rxn.GetReactantTemplate(i) for i in range(n_sub_t)]
    prod_mols = [rxn.GetProductTemplate(j) for j in range(n_prod_t)]

    # ATOM-UNBALANCED: a molecule was dropped before mapping (no SMILES for a
    # generic), so the mapper re-routed orphan atoms onto whatever remained. The
    # counts survive that; the pairs do not.
    if n_sub_t != len(sub_mnxms) or n_prod_t != len(prod_mnxms):
        return {}, "stripped"

    sub_named, amb_s = match_mols(sub_mols, sub_mnxms, canon)
    prod_named, amb_p = match_mols(prod_mols, prod_mnxms, canon)
    if amb_s or amb_p:
        return {}, "ambiguous_duplicate"

    sub_ranks = [canonical_ranks(m) for m in sub_mols]
    prod_ranks = [canonical_ranks(m) for m in prod_mols]
    if any(r is None for r in sub_ranks) or any(r is None for r in prod_ranks):
        return {}, "unrankable"

    # atom-map number -> (which substrate template, CANONICAL RANK, element)
    sub_index = {}
    for i, mol in enumerate(sub_mols):
        for a in mol.GetAtoms():
            n = a.GetAtomMapNum()
            if n > 0:
                sub_index[n] = (i, sub_ranks[i][a.GetIdx()], a.GetSymbol())

    pairs = defaultdict(list)
    for j, mol in enumerate(prod_mols):
        for a in mol.GetAtoms():
            n = a.GetAtomMapNum()
            if n <= 0:
                continue
            el = a.GetSymbol()
            if el not in ELEMENTS:
                continue
            src = sub_index.get(n)
            if src is None:
                continue
            si, s_rank, sel = src
            if sel != el:
                continue
            sm = sub_named[si] if si < len(sub_named) else None
            pm = prod_named[j] if j < len(prod_named) else None
            if not sm or not pm:
                continue
            # THE LINE THE UNIFIER DOES NOT WRITE: keep the pair, not two counts.
            pairs[(el, sm, pm)].append((s_rank, prod_ranks[j][a.GetIdx()]))
    if not pairs:
        return {}, "no_pairs"
    return dict(pairs), "ok"


# =====================================================================
# Driver
# =====================================================================

def cmd_extract(args):
    aam = pd.read_csv(args.aam, sep="\t")
    if args.min_confidence is not None and "confidence" in aam.columns:
        aam = aam[aam.confidence >= args.min_confidence]
    mnxrs = set(aam.mnxr)
    print(f"[atom-pairs] {len(aam):,} mapped reactions from {Path(args.aam).name}",
          flush=True)

    eqs = load_equations(Path(args.reac_prop), mnxrs)
    want = set()
    parsed = {}
    for r, eq in eqs.items():
        pe = parse_equation(eq)
        if pe:
            parsed[r] = pe
            want |= set(pe[0]) | set(pe[1])
    print(f"[atom-pairs] {len(parsed):,} equations parsed; {len(want):,} metabolites",
          flush=True)

    raw = load_mnxm_smiles(Path(args.chem_prop), want)
    canon = {}
    for m, smi in raw.items():
        cs = canon_smiles(smi)
        if cs:
            canon[m] = cs
    print(f"[atom-pairs] {len(canon):,} metabolites with canonical SMILES", flush=True)

    rows, status_rows = [], []
    tally = defaultdict(int)
    for k, rec in enumerate(aam.itertuples(index=False)):
        r = rec.mnxr
        pe = parsed.get(r)
        if pe is None:
            tally["no_equation"] += 1
            status_rows.append(dict(mnxr=r, status="no_equation", n_sub=0, n_prod=0,
                                    n_mapped_sub=0, n_mapped_prod=0))
            continue
        subs, prods = pe
        pairs, status = pairs_from_mapped(
            getattr(rec, "mapped_rxn_smiles", None), subs, prods, canon)
        tally[status] += 1
        status_rows.append(dict(mnxr=r, status=status, n_sub=len(subs),
                                n_prod=len(prods),
                                n_mapped_sub=len({p[1] for p in pairs}),
                                n_mapped_prod=len({p[2] for p in pairs})))
        for (el, sm, pm), idxs in pairs.items():
            rows.append(dict(mnxr=r, element=el, substrate=sm, product=pm,
                             n_atoms=len(idxs),
                             sub_idx=",".join(str(i) for i, _ in idxs),
                             prod_idx=",".join(str(j) for _, j in idxs)))
        if (k + 1) % 10000 == 0:
            print(f"[atom-pairs]   {k+1:,}/{len(aam):,}", flush=True)

    df = pd.DataFrame(rows, columns=list(PAIR_COLS))
    sdf = pd.DataFrame(status_rows, columns=list(STATUS_COLS))
    df.to_parquet(args.out, index=False)
    sdf.to_csv(args.out_status, sep="\t", index=False)

    print(f"\n[atom-pairs] {len(df):,} (rxn, element, sub, prod) pairs "
          f"over {df.mnxr.nunique():,} reactions -> {args.out}")
    print("\n[atom-pairs] per-reaction outcome:")
    for s, n in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"    {s:22s} {n:>7,}  ({100*n/len(aam):5.1f}%)")
    if len(df):
        print("\n[atom-pairs] pairs per element:")
        for el, g in df.groupby("element"):
            print(f"    {el}: {len(g):>7,} pairs  {g.n_atoms.sum():>9,} atoms  "
                  f"{g.mnxr.nunique():>6,} reactions")
    return 0


def cmd_selftest(args):
    """Does the pairing say what chemistry says? One reaction, checked by hand.

    MNXR106432: acetyl-CoA + CO2 + NADPH = pyruvate + NADP+ + CoA. The acetyl
    group's 2 carbons must land in pyruvate and the CoA moiety's 21 must land in
    CoA -- and pyruvate must receive nothing at all from NADPH.
    """
    smi = ("CC(=O)[S:1][CH2:2][CH2:3][NH:4][C:5](=[O:6])[CH2:7][CH2:8]"
           ">>[CH3:9][C:10](=[O:11])[C:12](=[O:13])[OH:14]")
    rxn = AllChem.ReactionFromSmarts(smi, useSmiles=True)
    print(f"parsed: {rxn is not None} "
          f"({rxn.GetNumReactantTemplates()} sub / {rxn.GetNumProductTemplates()} prod)")
    print("selftest is a parse check only; the real check is `verify` against MNX.")
    return 0


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("extract"); p.set_defaults(fn=cmd_extract)
    p.add_argument("--aam", required=True,
                   help="cached mapper output: mnxr, mapped_rxn_smiles, confidence")
    p.add_argument("--reac-prop", required=True)
    p.add_argument("--chem-prop", required=True)
    p.add_argument("--out", required=True, help="pairs parquet")
    p.add_argument("--out-status", required=True, help="per-reaction outcome tsv")
    p.add_argument("--min-confidence", type=float, default=None,
                   help="omit to keep the mapper's whole universe (tier B applies no "
                        "confidence gate; the acetyl-CoA split is chemically correct "
                        "at confidence 0.032, so a gate here would cost real pairs)")

    p = sub.add_parser("selftest"); p.set_defaults(fn=cmd_selftest)
    return ap.parse_args()


if __name__ == "__main__":
    a = parse_args()
    sys.exit(a.fn(a))
