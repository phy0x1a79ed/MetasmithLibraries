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

THE CROSSWALK THAT PARAGRAPH ASKS FOR: `--resolved`
---------------------------------------------------
The paragraph above names its own remedy and then declines to build it. `--resolved` is
that crosswalk: hand-authored, one row per claim, staged as a content-hashed input.

A RESOLVED metabolite is not a placeholder, and the difference is the whole design:

  * A PLACEHOLDER is scaffolding. It stands in for a carrier so the mapper sees a sane
    reaction, and its atoms are SUPPRESSED downstream -- they are not that molecule's
    atoms, so they must never become graph nodes.
  * A RESOLVED metabolite IS the metabolite. Its structure is supplied where MetaNetX has
    none, and its atoms are REAL: they become nodes, and they must, because the entire
    payoff is that the carrier's sulfur reaches the sink. Suppressing it would kill the
    very edge the row exists to license.

So a resolved structure joins the SMILES map and takes the ordinary path. Nothing special
happens to it downstream, which is deliberate: `ecspr_atom_pairs` needs no new filter,
because a metabolite that has a structure is what that module already handles.

WHAT MAKES IT HONEST, GIVEN THAT IT CANNOT BE VERIFIED
------------------------------------------------------
A curated row is an ASSERTION. It cannot be checked against MetaNetX -- "structureless"
means precisely that there is nothing there to check against (0/1531 blocking stubs carry
a formula, InChI, InChIKey or SMILES; 100% carry a name). Pretending otherwise would make
the name both the claim and its only support, which is trap #1 and has already produced
false verdicts in this lane. The warrant is the row's `basis` citation. The code's job is
not to validate the biology but to refuse the ways a row could be silently wrong:

  * ADDITIVE ONLY (`load_resolved`). A row may only supply a structure MetaNetX lacks;
    overriding an existing one is refused outright. That also means the crosswalk cannot
    perturb any reaction that maps today.
  * STALE-ID TRIPWIRE (`load_resolved`). The id is the key, and the recorded name must
    still match chem_prop. This is not name-as-proof -- it is the check that the id the
    curator reasoned about is the id they wrote down.
  * THE ROW CHECKS ITSELF (`load_resolved`). The asserted atom count must equal what the
    asserted SMILES actually contains.
  * BODIES MUST CANCEL (`gate_bodies_cancel`).

The `*` body is the honesty, not a shortcut: it says a carrier body exists and this row
does not claim to know it. Because the carrier appears on both sides (-SH in, -H out) the
unknown body cancels, so `concrete_balance` tests exactly the DIFFERENCE the row asserts --
one sulfur -- against the concrete chemistry either side of it. This is why B_verdict's
objection ("the body cancels, so a wrong carrier balances as well as the right one") kills
IDENTITY assertion but not this one: identity was never tested by balance, whereas the
difference is both what is claimed and what is tested. Where the bodies do NOT pair up the
cancellation argument fails, and `gate_bodies_cancel` refuses the reaction rather than
letting a `*`-counted-as-zero masquerade as a real count.
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
    (r"^\[?oxidized \[?adrenodoxin", FE_OX, "adrenodoxin"),
    (r"^\[?reduced \[?adrenodoxin", FE_RED, "adrenodoxin"),
    # diflavin NADPH--hemoprotein (cytochrome P450) reductase: a 2e- carrier via FAD/FMN,
    # no transferable C/N/S/P. Both bracket orderings occur ("oxidized [NADPH--...]" and
    # "[Oxidized NADPH---...]"), and MetaNetX writes the dash run as -- or ---.
    (r"^\[?oxidized \[?nadph[- ]+hemoprotein reductase", FE_OX, "hemoprotein_reductase"),
    (r"^\[?reduced \[?nadph[- ]+hemoprotein reductase", FE_RED, "hemoprotein_reductase"),
    # membrane cytochrome b5 (and other ferri/ferro-named cytochromes): ferri = oxidized (Fe3+),
    # ferro = reduced (Fe2+). One-electron protein carrier, no transferable C/N/S/P.
    (r"^ferricytochrome", FE_OX, "cytochrome_ferri"),
    (r"^ferrocytochrome", FE_RED, "cytochrome_ferri"),
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
    # generic redox donor/acceptor pair (A + 2[H] <-> AH2): structure unknown, transfers only
    # H, and a reaction naming it is a generic template, not a concrete instance.
    r"|a|ah2|reduced acceptor"
    # non-molecules: an electron and a photon carry no atom -- nothing for the atom lane.
    r"|e\(-\)|e-|hnu|hn|h\N{GREEK SMALL LETTER NU}|photon|light"
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

def _count_struct(smiles: str, X: str):
    """Atoms of element X in a SMILES, counted from the STRUCTURE.

    The `*` dummy contributes to no element -- which is the correct reading of a curated
    carrier, not a convenient one: the row asserts the drawn atoms and explicitly declines
    to claim the body. `gate_bodies_cancel` is what makes that silence safe for balance.
    """
    from rdkit import Chem, RDLogger      # imported here: diagnostics run this module
    RDLogger.DisableLog("rdApp.*")        # under envs without rdkit (p312)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() == X)


def load_resolved(path: Path, chem: dict):
    """mnxm -> curated SMILES. Every row is an ASSERTION; see the module docstring.

    The gates here cannot check the biology -- nothing can, which is why the row carries a
    citation. They refuse the ways a row could be silently wrong: overriding a structure
    MetaNetX already has, naming an id that no longer means what the curator thought, or
    asserting an atom count its own SMILES does not contain.
    """
    # `#` is a comment ONLY at line start. An inline `comment="#"` truncates any curated
    # SMILES that contains a `#` triple bond (nitrile C#N, alkyne C#C) -- silently dropping
    # exactly the rows a curator most needs to supply. Filter full-line comments and parse
    # the remainder, so a `#` inside a field survives.
    from io import StringIO
    _kept = [ln for ln in Path(path).read_text().splitlines() if not ln.lstrip().startswith("#")]
    d = pd.read_csv(StringIO("\n".join(_kept)), sep="\t")
    need = {"mnxm", "smiles", "mnx_name", "element", "n_atoms", "basis"}
    missing = need - set(d.columns)
    if missing:
        raise SystemExit(f"[resolved] crosswalk lacks columns: {sorted(missing)}")
    out = {}
    for r in d.itertuples(index=False):
        rec = chem.get(r.mnxm)
        if rec is None:
            raise SystemExit(f"[resolved] {r.mnxm}: not in chem_prop")
        nm, _formula, smi = rec
        if smi:
            raise SystemExit(
                f"[resolved] {r.mnxm} already HAS a structure ({smi!r}). A curated row may"
                f" only SUPPLY a structure MetaNetX lacks, never override one -- that is a"
                f" different and far larger claim.")
        if (nm or "").strip() != str(r.mnx_name).strip():
            raise SystemExit(
                f"[resolved] {r.mnxm}: row says {str(r.mnx_name)!r}, chem_prop says "
                f"{nm!r}. The id is the key -- a mismatch means the row is about a "
                f"different metabolite than the curator reasoned about.")
        got = _count_struct(r.smiles, r.element)
        if got is None:
            raise SystemExit(f"[resolved] {r.mnxm}: SMILES {r.smiles!r} does not parse")
        if got != int(r.n_atoms):
            raise SystemExit(
                f"[resolved] {r.mnxm}: asserts {r.n_atoms} {r.element}, but its SMILES "
                f"{r.smiles!r} contains {got}")
        # `pd.isna` FIRST: an empty cell arrives as NaN, and `str(nan)` is "nan" -- truthy,
        # so a bare `if not str(r.basis).strip()` accepts a row with no citation at all.
        # Measured: that gate silently passed the empty-basis case it was written to catch.
        if pd.isna(r.basis) or not str(r.basis).strip():
            raise SystemExit(f"[resolved] {r.mnxm}: no basis. A curated row cannot be "
                             f"verified against MetaNetX -- the citation IS its evidence.")
        out[r.mnxm] = r.smiles
    return out


def gate_bodies_cancel(subs, prods, resolved: dict):
    """Do the curated `*` bodies pair across the reaction?

    A curated carrier draws its body as `*` and counts it as zero for every element. That
    is only safe for balance when the same body stands on both sides and cancels. A carrier
    appearing on ONE side would have its unknown body silently counted as nothing, and the
    balance verdict would be about a molecule that does not exist.

    Counts `*` atoms among resolved participants per side and requires equality. It does
    NOT check the bodies are the same body -- it cannot; that is part of what the row
    asserts. It refuses the case where they demonstrably cannot cancel.
    """
    def n_star(ms):
        return sum(resolved[m].count("*") for m in ms if m in resolved)
    return n_star(subs) == n_star(prods)


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


def concrete_balance(subs, prods, chem, ph_mnxms, X, resolved=None):
    """Do the CONCRETE (non-placeholder) atoms of element X balance?

    This is what tests the conservation claim. If a carrier actually donated or
    absorbed an X atom, the concrete side that gained or lost it no longer balances,
    and this returns False -- so the reaction is refused for X while remaining usable
    for the elements the carrier really is inert to.

    Returns None when a concrete participant's formula is untrustworthy, which is a
    refusal too: an unknown count cannot be balanced.

    A RESOLVED participant is counted from its CURATED STRUCTURE, not skipped. It is not
    scaffolding -- it is a metabolite whose structure this run supplies, so its atoms are
    part of the chemistry being balanced, and its asserted difference is precisely what
    the balance then tests. Without `resolved` this behaves exactly as before: a curated
    carrier has no formula, `count_el` returns None, and the whole reaction would abstain
    -- the permissive answer, and the wrong one, since the gate is the point.
    """
    resolved = resolved or {}
    tot = {}
    for side, ms in (("s", subs), ("p", prods)):
        n = 0
        for m in ms:
            if m in ph_mnxms:
                continue
            if m in resolved:
                c = _count_struct(resolved[m], X)
            else:
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

    # Curated structures join the SMILES map, so a resolved metabolite stops being a
    # "generic" at `gens` below and simply becomes a metabolite with a structure. That is
    # the whole integration: no placeholder, no suppression, no downstream special case.
    resolved = load_resolved(Path(args.resolved), chem) if args.resolved else {}
    if resolved:
        print(f"\n[aam-rescue] {len(resolved)} CURATED structures -- each an ASSERTION, "
              f"warranted by its citation, not derived:")
        for m, s in resolved.items():
            print(f"    {m:<14} {str(chem[m][0])[:44]:<44} -> {s}")
        smi.update(resolved)

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
        # `smi` now carries the curated structures, so a reaction blocked ONLY by resolved
        # carriers has no generics left -- and would fall through here as "already
        # buildable" and never be mapped, rescuing exactly nothing. It is not already
        # buildable: the universe skipped it precisely because those participants had no
        # structure at the time. Supplying one is what makes it a rescue case.
        has_resolved = any(m in resolved for m in parts)
        if not gens and not has_resolved:
            tally["already buildable (not a placeholder case)"] += 1
            continue
        res = {m: placeholder_for(chem.get(m, ("", "", ""))[0]) for m in gens}
        if any(v is None for v in res.values()):
            tally["has a generic with no admissible placeholder"] += 1
            continue
        if not gate_bodies_cancel(subs, prods, resolved):
            # A curated `*` body counts as zero for every element. That is only sound when
            # the same body stands on both sides; here it demonstrably does not, so the
            # balance verdict would describe a molecule that does not exist.
            tally["REFUSED: curated bodies do not cancel"] += 1
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
            b = concrete_balance(subs, prods, chem, set(ph_smi), X, resolved)
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
            b = concrete_balance(subs, prods, chem, set(ph_smi), X, resolved)
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
    b.add_argument("--resolved", default=None,
                   help="curated structure crosswalk (mnxm, smiles, mnx_name, element, "
                        "n_atoms, basis). Each row ASSERTS a structure MetaNetX lacks, "
                        "warranted by its citation -- it is not derived and cannot be "
                        "auto-verified. Resolved metabolites are NOT placeholders: their "
                        "atoms are real and become graph nodes.")
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
