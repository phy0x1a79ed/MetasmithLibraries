"""Recover the reactions the AAM universe never mapped, by standing in for generic carriers.

WHAT THIS FIXES, AND WHY IT IS NOT THE EXTRACTOR'S PROBLEM
---------------------------------------------------------
`19b_rxnmapper_universe.py:build_rxn_smiles` needs a SMILES for EVERY participant:

    s = smi_map.get(m)
    if not s:
        return None, f"no_smiles:{m}"      # <- the whole reaction is discarded

So ONE structureless participant sinks the reaction before the mapper ever runs. That is
not a refusal tier `ecspr_atom_pairs.py` can see -- the reaction is simply absent from
the universe. Measured against the host (epi300) evidence set: 28.8% of host reactions
are missing this way, and 98% of those are missing for exactly this reason.

The casualties are not exotic. Sulfate assimilation is severed at sulfite -> H2S; the
MEP pathway is severed at HMBPP -> IPP, which is the ONLY reaction joining IPP to the
rest of carbon metabolism. In each case the structureless participant is an electron
carrier that contributes no atom to the chemistry at issue.

THE STAND-IN, AND THE CLAIM IT MAKES
------------------------------------
A placeholder is a deterministic, structure-bearing stand-in for a named generic. It
exists ONLY so the mapper sees a chemically sane reaction; its atoms are never emitted
as graph nodes (`--placeholders` makes the extractor drop any pair touching one).

Substituting one makes a real claim -- that the carrier is CONSERVED with respect to the
element being scored: every atom of element X it brings reappears in its partner, and
none crosses to a concrete metabolite. For ferredoxin in IspH that is true of carbon and
false of nothing. For the SAME ferredoxin in biotin synthase it is FALSE of sulfur: there
the [2Fe-2S] cluster IS the sulfur donor.

So the claim is NOT asserted per carrier -- it is TESTED per reaction and per element, by
`concrete_balance()`. A placeholder standing in for a carrier that actually donated an
atom leaves the concrete atoms of that element unbalanced, and the reaction is refused
for that element alone. Ferredoxin therefore rescues carbon in IspH and is correctly
refused for sulfur in BioB, with no special-casing anywhere.

WHAT DELIBERATELY GETS NO PLACEHOLDER
-------------------------------------
Three classes are refused outright, because a stand-in for them would be an invention
rather than a substitution:

  * Non-molecules. Participants named `Unknown`, `Carbon`, `d`, `IDH1`,
    `Enzyme-ligand complex`. There is no structure to stand in FOR. Note that
    MNXR180436 (`CMP + PG + 2 CoA + 2 phosphate = CTP + 2 "Carbon" + 2 glycerol-3-P`)
    is a SEED lump of five EC numbers -- 2.3.1.15, 2.3.1.51, 2.7.7.41, 2.7.8.5,
    3.1.3.27 -- and is the only host route to phosphatidylglycerol's carbon. Mapping it
    would manufacture the very edge the axis is meant to test.
  * Acyl carriers (ACP, holo-[ACP], malonyl-ACP, acyl-[ACP]). Unlike an electron
    carrier, an ACP's thioester DOES carry the atoms through. A stand-in would have to
    invent where they attach, and the balance gate cannot catch that -- the atoms would
    balance while being routed through fabricated bonds.
  * Unreconciled MNXref stubs -- `hexadecenoate` (seedM:cpd15237), `D-Isocitrate`
    (sabiorkM:30734), `Beta-D-Fructose 6-phosphate` (biggM:f6p_B). These have no
    formula, no InChI and no SMILES: MetaNetX minted an ID for a source compound it
    never reconciled to a structure. They are ordinary metabolites and the right fix is
    to resolve them to their structured entry, NOT to stand in for them. Resolving by
    name is what this codebase has already been burned by (a gate "failed" because
    MetaNetX calls it `CoA`, not `coenzyme A`), so they stay refused until there is an
    identifier-based crosswalk.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

EQ_TERM = re.compile(r"(\d+(?:\.\d+)?)\s+(MNXM\w+)@\w+")
SMILES_LEN_LIMIT = 8000
ELEMENTS = ("C", "N", "S", "P")

# ---------------------------------------------------------------------------
# The placeholder library. Order matters: first match wins.
#
# Each ox/red partner pair is atom-matched on purpose -- identical heavy-atom skeletons
# differing only in oxidation state -- so the mapper pairs them with each other trivially
# and cannot confuse them with a concrete metabolite.
# ---------------------------------------------------------------------------
FE_OX, FE_RED = "[Fe+3]", "[Fe+2]"
# thioredoxin-family: two cysteine thiols <-> one disulfide. Both C3S2; they differ by
# the 2 H that oxidation removes, which is the real chemistry.
DITHIOL, DISULFIDE = "SCCCS", "C1CCSS1"
# quinone pool: redox conserves every carbon of the ring. Both C6O2.
QUINONE, QUINOL = "O=C1C=CC(=O)C=C1", "Oc1ccc(O)cc1"

PLACEHOLDERS = [
    # --- one-electron / two-electron protein carriers: no transferable C, N, S or P ---
    (r"^oxidized \[?2fe-2s\]?[- ]\[?ferredoxin", FE_OX, "fe_s_carrier"),
    (r"^reduced \[?2fe-2s\]?[- ]\[?ferredoxin", FE_RED, "fe_s_carrier"),
    (r"^oxidized \[?4fe-4s\]?[- ]\[?ferredoxin", FE_OX, "fe_s_carrier"),
    (r"^reduced \[?4fe-4s\]?[- ]\[?ferredoxin", FE_RED, "fe_s_carrier"),
    (r"^oxidized \[?ferredoxin", FE_OX, "fe_s_carrier"),
    (r"^reduced \[?ferredoxin", FE_RED, "fe_s_carrier"),
    (r"^oxidized \[?flavodoxin", FE_OX, "flavo_carrier"),
    (r"^reduced \[?flavodoxin", FE_RED, "flavo_carrier"),
    (r"^oxidized \[?electron-transfer flavoprotein", FE_OX, "etf_carrier"),
    (r"^reduced \[?electron-transfer flavoprotein", FE_RED, "etf_carrier"),
    (r"^fe\(iii\)-\[?cytochrome", FE_OX, "cytochrome"),
    (r"^fe\(ii\)-\[?cytochrome", FE_RED, "cytochrome"),
    (r"^oxidized \[?rubredoxin", FE_OX, "rubredoxin"),
    (r"^reduced \[?rubredoxin", FE_RED, "rubredoxin"),
    # --- thiol/disulfide redox carriers: sulfur-bearing, and sulfur-conserving ---
    (r"^\[?thioredoxin\]?-dithiol", DITHIOL, "thiol_carrier"),
    (r"^\[?thioredoxin\]?-disulfide", DISULFIDE, "thiol_carrier"),
    (r"^\[?glutaredoxin\]?-dithiol", DITHIOL, "thiol_carrier"),
    (r"^\[?glutaredoxin\]?-disulfide", DISULFIDE, "thiol_carrier"),
    # --- quinone pool ---
    (r"^an? ubiquinone$", QUINONE, "quinone"),
    (r"^an? ubiquinol$", QUINOL, "quinone"),
    (r"^an? menaquinone$", QUINONE, "quinone"),
    (r"^an? menaquinol$", QUINOL, "quinone"),
    (r"^an? demethylmenaquinone$", QUINONE, "quinone"),
    (r"^an? demethylmenaquinol$", QUINOL, "quinone"),
]
_COMPILED = [(re.compile(p, re.I), s, t) for p, s, t in PLACEHOLDERS]

# Never stand in for these, whatever else matches. See the module docstring.
REFUSE = re.compile(
    r"^(unknown|carbon|d|nad|nadh|idh\d*|enzyme-\w+ complex|acceptor|reduced acceptor"
    r"|.*\bacp\b.*|.*acyl-carrier.*|.*acyl carrier.*|starch|chitin|.*tRNA.*"
    r"|phosphoprotein|.*\[protein\]|protein .*)$", re.I)


def placeholder_for(name: str):
    """(smiles, tag) for a structureless generic, or None to leave it refused."""
    if not name:
        return None
    n = name.strip()
    if REFUSE.match(n):
        return None
    for rx, smi, tag in _COMPILED:
        if rx.match(n):
            return smi, tag
    return None


# ---------------------------------------------------------------------------
# MetaNetX inputs (same derivations the universe builder uses)
# ---------------------------------------------------------------------------

def load_chem(chem_prop: Path):
    """mnxm -> (name, formula, smiles)."""
    out = {}
    with open(chem_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 9:
                out[p[0]] = (p[1], p[3], p[8].strip())
    return out


def load_equations(reac_prop: Path, want=None):
    out = {}
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) > 1 and p[0] != "EMPTY" and (want is None or p[0] in want):
                out[p[0]] = p[1]
    return out


def parse_equation(eq: str):
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


def build_rxn_smiles(subs, prods, smi_map):
    """Identical in convention to 19b's builder: one fragment per EQ_TERM match, in
    equation order, joined by '.'. WATER/BIOMASS are not MNXM tokens so they get no
    fragment -- the extractor's regex must agree with that, and does."""
    def side(ms):
        out = []
        for m in ms:
            s = smi_map.get(m)
            if not s:
                return None, m
            out.append(s)
        return out, None

    Ls, bad = side(subs)
    if bad:
        return None, bad
    Rs, bad = side(prods)
    if bad:
        return None, bad
    return ".".join(Ls) + ">>" + ".".join(Rs), None


# ---------------------------------------------------------------------------
# The gate that makes the substitution honest
# ---------------------------------------------------------------------------
_FORM = re.compile(r"([A-Z][a-z]?)(\d*)")


def count_el(formula: str, X: str):
    """Atoms of element X in a MetaNetX formula; None when it cannot be trusted
    (absent, wildcard `*` polymer, or nested groups)."""
    if not formula or formula.strip() in ("", "*"):
        return None
    if "*" in formula or "(" in formula or ")" in formula:
        return None
    n, seen = 0, False
    for sym, num in _FORM.findall(formula):
        if not sym:
            continue
        if sym == X:
            seen = True
            n += int(num) if num else 1
    return n if seen else 0


def concrete_balance(subs, prods, chem, ph_mnxms, X):
    """Do the CONCRETE (non-placeholder) atoms of element X balance?

    This is what tests the conservation claim. If a carrier actually donated or
    absorbed an X atom, the concrete side that gained or lost it no longer balances,
    and this returns False -- so the reaction is refused for X while remaining usable
    for the elements the carrier really is inert to.

    Returns None when a concrete participant's formula is untrustworthy, which is a
    refusal too: an unknown count cannot be balanced.
    """
    tot = {}
    for side, ms in (("s", subs), ("p", prods)):
        n = 0
        for m in ms:
            if m in ph_mnxms:
                continue
            c = count_el(chem.get(m, ("", "", ""))[1], X)
            if c is None:
                return None
            n += c
        tot[side] = n
    if tot["s"] == 0 and tot["p"] == 0:
        return None          # element absent; nothing to say
    return tot["s"] == tot["p"]


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def cmd_build(args):
    chem = load_chem(Path(args.chem_prop))
    smi = {m: v[2] for m, v in chem.items() if v[2]}
    print(f"[aam-rescue] chem_prop: {len(chem):,} metabolites, {len(smi):,} with SMILES",
          flush=True)

    want = None
    if args.targets:
        want = {l.strip() for l in open(args.targets) if l.strip()}
        print(f"[aam-rescue] restricted to {len(want):,} target reactions", flush=True)
    skip = set()
    if args.exclude:
        d = pd.read_csv(args.exclude, sep="\t")
        skip = set(d[d.columns[0]])
        print(f"[aam-rescue] excluding {len(skip):,} already-mapped reactions", flush=True)

    eqs = load_equations(Path(args.reac_prop), want)
    print(f"[aam-rescue] {len(eqs):,} equations", flush=True)

    ph_smi, ph_tag = {}, {}
    todo, tally = [], Counter()
    for r, eq in eqs.items():
        if r in skip:
            continue
        pe = parse_equation(eq)
        if not pe:
            tally["unparseable equation"] += 1
            continue
        subs, prods = pe
        parts = set(subs) | set(prods)
        gens = [m for m in parts if m not in smi]
        if not gens:
            tally["already buildable (not a placeholder case)"] += 1
            continue
        res = {m: placeholder_for(chem.get(m, ("", "", ""))[0]) for m in gens}
        if any(v is None for v in res.values()):
            tally["has a generic with no admissible placeholder"] += 1
            continue
        for m, (s, t) in res.items():
            ph_smi[m], ph_tag[m] = s, t
        todo.append((r, subs, prods))
        tally["PLACEHOLDER-RESOLVABLE"] += 1

    print("\n[aam-rescue] triage:")
    for k, n in tally.most_common():
        print(f"    {k:<48} {n:>6,}")
    print(f"\n[aam-rescue] {len(ph_smi):,} distinct generics stood in for:")
    for m, s in sorted(ph_smi.items(), key=lambda x: ph_tag[x[0]]):
        print(f"    {m:<14} {str(chem[m][0])[:44]:<44} [{ph_tag[m]:<14}] -> {s}")

    merged = dict(smi)
    merged.update(ph_smi)
    rows = []
    for r, subs, prods in todo:
        s, bad = build_rxn_smiles(subs, prods, merged)
        if not s:
            continue
        if len(s) > SMILES_LEN_LIMIT:
            continue
        rows.append((r, s, subs, prods))
    print(f"\n[aam-rescue] {len(rows):,} reactions built into SMILES", flush=True)
    if not rows:
        print("[aam-rescue] nothing to map")
        return

    from rxnmapper import RXNMapper
    rm = RXNMapper()
    out, bs = [], args.batch_size
    for i in range(0, len(rows), bs):
        chunk = rows[i:i + bs]
        try:
            res = rm.get_attention_guided_atom_maps([c[1] for c in chunk])
        except Exception as e:
            print(f"[aam-rescue]   batch {i} failed ({e}); falling back to singles",
                  flush=True)
            res = []
            for c in chunk:
                try:
                    res.append(rm.get_attention_guided_atom_maps([c[1]])[0])
                except Exception:
                    res.append({"mapped_rxn": None, "confidence": 0.0})
        for c, rr in zip(chunk, res):
            out.append(dict(mnxr=c[0], rxn_smiles=c[1],
                            mapped_rxn_smiles=rr.get("mapped_rxn"),
                            confidence=rr.get("confidence", 0.0)))
        if (i // bs) % 20 == 0:
            print(f"[aam-rescue]   mapped {min(i+bs, len(rows)):,}/{len(rows):,}",
                  flush=True)

    df = pd.DataFrame(out, columns=["mnxr", "rxn_smiles", "mapped_rxn_smiles",
                                    "confidence"])
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\n[aam-rescue] {len(df):,} newly mapped reactions -> {args.out}")

    # the placeholder map: SMILES so the extractor can NAME these fragments, and the
    # per-element admissibility the balance gate just decided.
    prows = []
    eq_by_r = {r: (s, p) for r, _, s, p in rows}
    ok_el = {X: 0 for X in ELEMENTS}
    for r, subs, prods in [(r, s, p) for r, _, s, p in rows]:
        for X in ELEMENTS:
            b = concrete_balance(subs, prods, chem, set(ph_smi), X)
            if b:
                ok_el[X] += 1
    print("[aam-rescue] reactions whose CONCRETE atoms balance, per element:")
    for X in ELEMENTS:
        print(f"    {X}: {ok_el[X]:,}/{len(rows):,}")

    for m, s in ph_smi.items():
        prows.append(dict(mnxm=m, smiles=s, tag=ph_tag[m], name=chem[m][0]))
    pd.DataFrame(prows).to_csv(args.out_placeholders, sep="\t", index=False)
    print(f"[aam-rescue] placeholder map -> {args.out_placeholders}")

    # per-reaction, per-element admissibility -- the extractor applies this as a filter
    brows = []
    for r, subs, prods in [(r, s, p) for r, _, s, p in rows]:
        for X in ELEMENTS:
            b = concrete_balance(subs, prods, chem, set(ph_smi), X)
            if b is not None:
                brows.append(dict(mnxr=r, element=X, balanced=bool(b)))
    pd.DataFrame(brows).to_csv(args.out_balance, sep="\t", index=False)
    print(f"[aam-rescue] per-element balance -> {args.out_balance}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="map reactions the universe skipped")
    b.add_argument("--reac-prop", required=True)
    b.add_argument("--chem-prop", required=True)
    b.add_argument("--targets", help="file of MNXR ids, one per line (default: all)")
    b.add_argument("--exclude", help="TSV whose first column lists already-mapped MNXRs")
    b.add_argument("--out", required=True)
    b.add_argument("--out-placeholders", required=True)
    b.add_argument("--out-balance", required=True)
    b.add_argument("--batch-size", type=int, default=8)
    b.set_defaults(func=cmd_build)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
