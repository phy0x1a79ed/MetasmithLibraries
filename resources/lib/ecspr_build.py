"""Build an :class:`ecspr_graph.AtomGraph` -- from atom pairs, from a GEM, or from a GPR.

Three builders, one output. Everything downstream measures the same object, so what a
number means is decided entirely by which reaction set went in and with what conductance.

  * :func:`graph_from_pairs` -- the primitive. Atom-transfer table x per-reaction evidence
    weight x per-reaction direction ratio -> a directed atom network. The other two reduce
    to it.
  * :func:`gem_to_graph`     -- a curated genome-scale model (COBRA JSON **or** SBML),
    crosswalked to canonical MNXR and induced on the atom-mapped universe with uniform
    E = 1.0.
  * :func:`gpr_to_graph`     -- a GEM plus a set of ACTIVE GENES. Reactions whose
    gene-protein-reaction rule is unsatisfied are dropped.

WHY THE GPR BUILDER IS THE POINT
--------------------------------
It is what lets the engine stop owning perturbation. A knockout is not an argument to the
solver; it is a smaller gene set, hence a different network, and the comparison is a
subtraction the caller does on two independent solves. This replaces the cross-scope
`_val.ko_dead_mnxrs` the method map documents, and it means the engine never learns what a
knockout is.

The GPR is evaluated by **cobra's own** ``Reaction.gpr.eval`` -- the only correct AND/OR
resolver. Re-deriving it from the gene-reaction rule string is how you silently turn an
isozyme pair into a required complex.

AAM COVERAGE IS A NUMBER, NOT A SILENCE
---------------------------------------
Reactions absent from the atom-pair table cannot contribute an edge. That set is the AAM
gap, and every builder returns its size in ``graph.meta``. It is the honest coverage of the
built network, and a build that does not report it is a build whose conductances cannot be
interpreted.

Env: numpy + pandas (+ cobra only on the GEM/GPR paths, imported lazily).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ecspr_graph import AtomGraph, Terminal, solve                     # noqa: F401

ELEMENTS = ("C", "N", "S", "P")


# =====================================================================
# Loaders
# =====================================================================

def load_pairs(path, element=None) -> pd.DataFrame:
    """The atom-transfer table. ``element`` restricts up front -- the carbon slice alone is
    ~2M rows, so filtering before anything else is the difference between seconds and
    minutes."""
    df = pd.read_parquet(path)
    return df if element is None else df[df.element == element]


def load_direction_ratios(path) -> dict:
    """``{mnxr: g_rev/g_fwd}`` from the direction ensemble. A reaction absent from the table
    has ratio 1.0 -- the symmetric limit, i.e. "no directional evidence" degrades exactly to
    the undirected model rather than to a guess."""
    df = pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path, sep="\t")
    return {str(r): float(v) for r, v in zip(df.mnxr, df.ratio) if str(r) != "EMPTY"}


def load_evidence_weights(path, source="epi300", column="E_full") -> dict:
    """``{mnxr: E_r}`` for one evidence source. ``E_full`` is the belief-conserving
    allocation; ``E_dlec`` is its domain-limited variant."""
    df = pd.read_parquet(path)
    if "source" in df.columns and source is not None:
        df = df[df.source == source]
    return {str(r): float(e) for r, e in zip(df.mnxr, df[column]) if float(e) > 0}


# =====================================================================
# The primitive builder
# =====================================================================

def _explode_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise the two atom-pair schemas onto one row per atom correspondence.

    The frozen reference is already at pair granularity (int ``sub_idx``/``prod_idx``, a
    scalar ``pair_w``). The incumbent extract and the toy fixtures pack a whole
    substrate->product atom list into one row as comma-joined strings. Both must keep
    working; the same two-schema tolerance ``ecspr_atom_graph._rank_list`` /
    ``_weight_list`` carry, expressed once as a reshape instead of per-row in a loop.
    """
    if not (df.sub_idx.dtype == object or df.prod_idx.dtype == object):
        return df
    rows = []
    for rec in df.itertuples(index=False):
        si = [int(x) for x in str(rec.sub_idx).split(",")]
        pi = [int(x) for x in str(rec.prod_idx).split(",")]
        raw = getattr(rec, "pair_w", None)
        if isinstance(raw, str) and raw:
            pw = [float(x) for x in raw.split(",")]
        else:
            try:
                f = float(raw)
                pw = [1.0 if f != f else f] * len(si)
            except (TypeError, ValueError):
                pw = [1.0] * len(si)
        for a, b, w in zip(si, pi, pw):
            rows.append((rec.mnxr, rec.element, rec.substrate, rec.product, a, b, w,
                         getattr(rec, "confidence", 1.0)))
    return pd.DataFrame(rows, columns=["mnxr", "element", "substrate", "product",
                                       "sub_idx", "prod_idx", "pair_w", "confidence"])


def graph_from_pairs(pairs: pd.DataFrame, element: str, weights: dict,
                     ratios: dict | None = None, *, use_confidence: bool = False,
                     with_provenance: bool = False,
                     meta: dict | None = None) -> AtomGraph:
    """Directed atom network for one element.

    An edge is one atom transfer: tail ``(substrate, sub_idx)``, head ``(product,
    prod_idx)``, with ``gp = E_r * pair_w`` and ``gm = ratio_r * gp``. Orientation comes
    straight from the pair table's substrate/product columns -- no reaction roles, no
    reac_prop parsing, because on an atom graph the transfer already names its own
    direction.

    ``pair_w`` is the pair's fanout-dilution weight: 1.0 for an ensemble consensus, 0.5 for
    a lone-member correspondence, and conf/wsum for a disagreement diluted across
    candidates. It is margin-conserving (a source atom's candidate weights sum to 1.0),
    which is why it -- and not the raw ``confidence`` column -- is what folds into the
    conductance. ``use_confidence`` multiplies ``confidence`` in as a SENSITIVITY knob; it
    double-counts the dilution on lone-member and disagreement rows, so it is off by
    default.

    Parallel edges in the SAME direction are summed. That is exact, not an approximation:
    both see the same potential drop ``x``, so ``sum_i [gp_i h(x) + gm_i (x - h(x))]`` is
    the single edge ``(sum gp_i, sum gm_i)``. Anti-parallel edges are kept separate, where
    the same identity does not hold.

    A RATIO ABOVE 1 FLIPS THE EDGE -- it does not amplify it
    -------------------------------------------------------
    ``ratio = g_rev/g_fwd = exp(dG'/RT)`` is defined against the direction the MetaNetX
    equation is WRITTEN in, so ``ratio > 1`` means the reaction as written is
    thermodynamically unfavourable and actually runs the other way. On the real ensemble
    that is not a rare edge case: ratios span 3e-18 to 3.3e17, and a quarter of the table
    sits above 1.

    Taking ``gm = ratio * gp`` literally there would give the reverse branch a conductance
    of ``3e17 * E_r``, i.e. it would let the DIRECTION evidence manufacture conductance the
    EXISTENCE evidence never supported -- a short circuit, and the reason a first run of
    this builder returned per-precursor currents ~1e13 times the injected current. Direction
    evidence may only ever throttle, never amplify. So an edge with ``ratio > 1`` is
    REVERSED (tail and head swapped) and its ratio inverted, which preserves the physical
    forward/backward asymmetry exactly while keeping the favoured branch at ``E_r *
    pair_w``. At ``ratio == 1`` it is the identity, so the symmetric limit is untouched.
    """
    ratios = ratios or {}
    df = pairs[pairs.element == element] if "element" in pairs.columns else pairs
    df = df[df.mnxr.isin(weights.keys())]
    df = _explode_schema(df)
    n_rxn_in = len(weights)
    if df.empty:
        return AtomGraph([], [], np.zeros(0), np.zeros(0),
                         dict(meta or {}, element=element, n_reactions_requested=n_rxn_in,
                              n_reactions_used=0, n_aam_gap=n_rxn_in, n_nodes=0, n_edges=0))

    er = df.mnxr.map(weights).astype(float).to_numpy()
    pw = pd.to_numeric(df.get("pair_w", 1.0), errors="coerce").fillna(1.0).to_numpy() \
        if "pair_w" in df.columns else np.ones(len(df))
    gp = er * pw
    if use_confidence:
        cf = pd.to_numeric(df.get("confidence", 1.0), errors="coerce").fillna(1.0).to_numpy() \
            if "confidence" in df.columns else np.ones(len(df))
        gp = gp * cf
    ratio = df.mnxr.map(ratios).astype(float).fillna(1.0).to_numpy()

    # ratio > 1 => the reaction runs against the way its equation is written: flip the edge
    # and invert the ratio, so the favoured branch keeps conductance E_r * pair_w and the
    # direction evidence can only ever throttle. See the docstring.
    flip = ratio > 1.0
    t_met = np.where(flip, df["product"].to_numpy(), df.substrate.to_numpy())
    t_idx = np.where(flip, df.prod_idx.to_numpy(), df.sub_idx.to_numpy())
    h_met = np.where(flip, df.substrate.to_numpy(), df["product"].to_numpy())
    h_idx = np.where(flip, df.sub_idx.to_numpy(), df.prod_idx.to_numpy())
    ratio = np.where(flip, 1.0 / np.maximum(ratio, np.finfo(float).tiny), ratio)

    gm_row = ratio * gp
    keep = (gp > 0) & ~((t_met == h_met) & (t_idx == h_idx))
    if not keep.any():
        return AtomGraph([], [], np.zeros(0), np.zeros(0),
                         dict(meta or {}, element=element, n_reactions_requested=n_rxn_in,
                              n_reactions_used=0, n_aam_gap=n_rxn_in, n_nodes=0, n_edges=0))
    t_met, t_idx, h_met, h_idx = t_met[keep], t_idx[keep], h_met[keep], h_idx[keep]
    gp, gm_row = gp[keep], gm_row[keep]
    mnxr_row = df.mnxr.to_numpy()[keep]

    # One edge per distinct ORDERED (tail, head) atom pair; parallel rows sum onto it.
    key = pd.MultiIndex.from_arrays([t_met, t_idx, h_met, h_idx])
    codes, uniq = pd.factorize(key, sort=False)
    ne = len(uniq)
    gp_e = np.bincount(codes, gp, minlength=ne)
    gm_e = np.bincount(codes, gm_row, minlength=ne)

    nodes, idx, edges = [], {}, []

    def _i(k):
        j = idx.get(k)
        if j is None:
            j = idx[k] = len(nodes)
            nodes.append(k)
        return j

    for tm, ti, hm, hi in uniq:
        edges.append((_i((tm, int(ti))), _i((hm, int(hi)))))

    used = set(df.mnxr.unique())
    m = dict(meta or {})
    m.update(element=element,
             n_reactions_requested=n_rxn_in,
             n_reactions_used=len(used),
             n_aam_gap=n_rxn_in - len(used),
             n_pair_rows=int(len(df)),
             n_nodes=len(nodes), n_edges=len(edges),
             n_metabolites=len({k[0] for k in nodes}),
             n_directed_rows=int((ratio != 1.0).sum()),
             n_reversed_rows=int(flip.sum()))
    if with_provenance:
        # Which reactions built which edge, and with how much of its conductance. Current
        # divides between parallel conductances in exact proportion to conductance, so this
        # is what makes per-reaction current attribution exact rather than a heuristic.
        m["edge_reactions"] = pd.DataFrame(dict(edge=codes, mnxr=mnxr_row, gp=gp))
    return AtomGraph(nodes, edges, gp_e, gm_e, m, idx)


def reaction_currents(graph: AtomGraph, solution) -> pd.Series:
    """Per-reaction current, ``{mnxr: sum_e |i_e| * gp_r,e / gp_e}``, descending.

    Requires a graph built with ``with_provenance=True``. Current through parallel
    conductances divides in exact proportion to conductance, so splitting an edge's current
    by each contributing reaction's share of ``gp`` is exact, not an approximation.

    This is how a validation run picks its knockout FROM THE RUN'S OWN RANKING rather than
    naming a reaction in advance -- a demo that names its reaction proves nothing once the
    universe moves.
    """
    prov = graph.meta.get("edge_reactions")
    if prov is None:
        raise ValueError("graph was not built with with_provenance=True")
    cur = np.zeros(graph.m)
    oe, ie = solution.edge_currents()
    np.add.at(cur, oe, np.abs(ie))
    pe = prov.edge.to_numpy()
    denom = graph.gp[pe]
    frac = prov.gp.to_numpy() / np.where(denom > 0, denom, 1.0)
    return (pd.Series(cur[pe] * frac, index=prov.mnxr.to_numpy())
            .groupby(level=0).sum().sort_values(ascending=False))


# =====================================================================
# GEM -> graph
# =====================================================================

def load_model(path):
    """A COBRA model from JSON, SBML (.xml/.sbml/.xml.gz) or MATLAB. ``import cobra`` stays
    lazy so the pairs-only path never pays for it."""
    import cobra
    p = Path(path)
    name = p.name.lower()
    if name.endswith(".json"):
        return cobra.io.load_json_model(str(p))
    if name.endswith((".xml", ".sbml", ".xml.gz", ".sbml.gz")):
        return cobra.io.read_sbml_model(str(p))
    if name.endswith(".mat"):
        return cobra.io.load_matlab_model(str(p))
    if name.endswith((".yml", ".yaml")):
        return cobra.io.load_yaml_model(str(p))
    raise ValueError(f"unrecognised model format: {p.name}")


def load_reac_xref(path) -> tuple:
    """``reac_xref.tsv`` -> ``(bigg.reaction -> MNXR, kegg.reaction -> MNXR)``, first wins."""
    bigg2m, kegg2m = {}, {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1].startswith("MNXR"):
                continue
            if p[0].startswith("bigg.reaction:"):
                bigg2m.setdefault(p[0].split(":", 1)[1], p[1])
            elif p[0].startswith("kegg.reaction:"):
                kegg2m.setdefault(p[0].split(":", 1)[1], p[1])
    return bigg2m, kegg2m


def atom_universe(pairs: pd.DataFrame, elements=ELEMENTS) -> set:
    """Every MNXR carrying at least one atom transfer on one of ``elements``.

    THIS is what "in universe" means in the atom lane. The star lane asked the same
    question of ``mnx_bipartite`` -- whether a reaction carried a ``w_X > 0`` edge -- which
    is a property of a topology that no longer exists here. Same question, asked of the
    table the edges actually come from.
    """
    df = pairs[pairs.element.isin(list(elements))] if "element" in pairs.columns else pairs
    return set(df.mnxr.unique())


def crosswalk_gem(model, reac_xref, universe: set) -> tuple:
    """GEM reaction -> ONE canonical current MNXR.

    Preference order: a candidate that is in the atom-mapped universe first, then by source
    priority ``bigg > kegg > embedded``. Returns ``(DataFrame[rxn_id, mnxr, source,
    in_universe], stats)``. Lifted from the star lane's ``gem_crosswalk`` with the universe
    supplied rather than read off a bipartite pickle.
    """
    bigg2m, kegg2m = load_reac_xref(reac_xref)

    def as_list(v):
        return [] if not v else (v if isinstance(v, list) else [v])

    rows, n_in, n_gap, n_unresolved = [], 0, 0, 0
    for r in model.reactions:
        ann = r.annotation or {}
        cands = []
        if r.id in bigg2m:
            cands.append((0, "bigg", bigg2m[r.id]))
        for kk in as_list(ann.get("kegg.reaction")):
            if kk in kegg2m:
                cands.append((1, "kegg", kegg2m[kk]))
        for xx in as_list(ann.get("metanetx.reaction")):
            cands.append((2, "embedded", xx))
        if not cands:
            n_unresolved += 1
            continue
        inu = [c for c in cands if c[2] in universe]
        pick = min(inu or cands, key=lambda c: c[0])
        ok = pick[2] in universe
        rows.append(dict(rxn_id=r.id, mnxr=pick[2], source=pick[1], in_universe=ok))
        n_in += ok
        n_gap += (not ok)
    df = pd.DataFrame(rows, columns=["rxn_id", "mnxr", "source", "in_universe"])
    if not df.empty:
        df = df.drop_duplicates("rxn_id")
    stats = dict(n_reactions=len(model.reactions), n_resolved=len(df),
                 n_unresolved=n_unresolved, n_in_universe=n_in, n_aam_gap=n_gap)
    return df, stats


def gem_to_graph(model_path, reac_xref, pairs: pd.DataFrame, element: str,
                 ratios: dict | None = None, *, E: float = 1.0) -> AtomGraph:
    """A curated GEM's reactome, induced on the atom-mapped universe with uniform ``E``.

    Uniform conductance is deliberate: a curated model asserts that a reaction is PRESENT,
    not how much evidence there is for it, so weighting it by anything would be inventing a
    quantity. The evidence-weighted lane is :func:`graph_from_pairs` with real ``E_r``.
    """
    model = load_model(model_path)
    xw, stats = crosswalk_gem(model, reac_xref, atom_universe(pairs))
    weights = {r: float(E) for r in xw.mnxr.dropna().unique()}
    meta = dict(builder="gem", model=str(Path(model_path).name), uniform_E=float(E))
    meta.update({f"crosswalk_{k}": v for k, v in stats.items()})
    g = graph_from_pairs(pairs, element, weights, ratios, meta=meta)
    g.meta["crosswalk"] = xw
    return g


# =====================================================================
# GPR + active gene set -> graph
# =====================================================================

def gpr_to_graph(model_path, genes_active, reac_xref, pairs: pd.DataFrame, element: str,
                 ratios: dict | None = None, *, E: float = 1.0,
                 keep_ruleless: bool = True) -> AtomGraph:
    """The GEM restricted to reactions whose GPR is satisfied by ``genes_active``.

    An MNXR is LIVE if **any** GEM reaction mapping to it is live. The crosswalk is
    many-GEM-reactions-to-one-MNXR, so a live and a dead reaction can land on the same
    MNXR; treating that MNXR as dead would delete a route the organism still has. Recorded
    here rather than left to the reader.

    ``keep_ruleless`` keeps reactions with an empty gene-reaction rule -- exchanges,
    diffusion, spontaneous chemistry. They have no gene to knock out, so dropping them
    would make every gene set look like a starvation.
    """
    model = load_model(model_path)
    active = set(genes_active)
    # cobra's `GPR.eval` takes KNOCKOUTS, not active genes -- it answers "is this rule
    # still satisfied after removing these?". Handing it the active set inverts the
    # question: a first run of this builder reported 532 of 2742 reactions live for the
    # WILD TYPE and MORE live after a knockout. The active set is the caller's contract
    # because that is what a condition is; the complement is computed here, once.
    knockouts = {g.id for g in model.genes} - active
    live_ids, dead_ids = set(), set()
    for r in model.reactions:
        rule = (r.gene_reaction_rule or "").strip()
        if not rule:
            (live_ids if keep_ruleless else dead_ids).add(r.id)
            continue
        try:
            ok = bool(r.gpr.eval(knockouts))
        except Exception:                       # malformed rule -> treat as ruleless
            ok = keep_ruleless
        (live_ids if ok else dead_ids).add(r.id)

    xw, stats = crosswalk_gem(model, reac_xref, atom_universe(pairs))
    live_mnxr = set(xw[xw.rxn_id.isin(live_ids)].mnxr.dropna().unique())
    all_mnxr = set(xw.mnxr.dropna().unique())
    weights = {r: float(E) for r in live_mnxr}

    meta = dict(builder="gpr", model=str(Path(model_path).name), uniform_E=float(E),
                n_genes_active=len(active), n_genes_model=len(model.genes),
                n_genes_knocked_out=len(knockouts),
                n_rxn_live=len(live_ids), n_rxn_dead=len(dead_ids),
                n_mnxr_live=len(live_mnxr), n_mnxr_dead=len(all_mnxr - live_mnxr))
    meta.update({f"crosswalk_{k}": v for k, v in stats.items()})
    g = graph_from_pairs(pairs, element, weights, ratios, meta=meta)
    g.meta["crosswalk"] = xw
    return g


def model_gene_ids(model_path) -> list:
    """Every gene id in a model -- the full active set, i.e. the wild type."""
    return [g.id for g in load_model(model_path).genes]


# =====================================================================
# Self-tests
# =====================================================================

GRAPH_STEM = "atom_graph_{X}.npz"


def write_graph_dir(out_dir, graphs: dict, extra: dict | None = None):
    """An ``ecspr::atom_graph`` directory: one ``atom_graph_{X}.npz`` per element plus a
    build manifest carrying reaction counts, the AAM gap and the element."""
    import json
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = dict(extra or {})
    for X, g in graphs.items():
        g.save(out_dir / GRAPH_STEM.format(X=X))
        manifest[X] = {k: v for k, v in g.meta.items()
                       if not hasattr(v, "shape") and k != "edge_reactions"}
    (out_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    return out_dir


def read_graph_dir(graph_dir, element):
    return AtomGraph.load(Path(graph_dir) / GRAPH_STEM.format(X=element))


def _elements_in(graph_dir) -> list:
    import re
    pat = re.compile(GRAPH_STEM.format(X=r"(?P<X>\w+)").replace(".", r"\.")
                     .replace(r"\.npz", r"\.npz$"))
    return sorted(m.group("X") for f in Path(graph_dir).iterdir()
                  if (m := pat.match(f.name)))


# =====================================================================
# CLI
# =====================================================================

def cmd_build_graph(args):
    ratios = load_direction_ratios(args.direction) if args.direction else {}
    pairs_all = load_pairs(args.atom_pairs)
    elements = args.elements or ELEMENTS
    graphs, extra = {}, dict(source=args.source)
    if args.source == "pairs":
        weights = load_evidence_weights(args.weights, source=args.evidence_source,
                                        column=args.evidence_column)
        extra.update(evidence_source=args.evidence_source,
                     evidence_column=args.evidence_column)
        for X in elements:
            graphs[X] = graph_from_pairs(pairs_all, X, weights, ratios)
    elif args.source == "gem":
        for X in elements:
            graphs[X] = gem_to_graph(args.model, args.reac_xref, pairs_all, X, ratios,
                                     E=args.uniform_e)
    else:                                                       # gpr
        genes = [ln.strip() for ln in Path(args.gene_set).read_text().split()
                 if ln.strip()]
        extra.update(n_genes_active=len(genes))
        for X in elements:
            graphs[X] = gpr_to_graph(args.model, genes, args.reac_xref, pairs_all, X,
                                     ratios, E=args.uniform_e)
    for X, g in graphs.items():
        print(f"[build-graph] {X}: {g.n:,} nodes / {g.m:,} edges from "
              f"{g.meta['n_reactions_used']:,} reactions (AAM gap "
              f"{g.meta['n_aam_gap']:,})", flush=True)
    write_graph_dir(args.out_dir, graphs, extra)
    print(f"[build-graph] -> {args.out_dir}")
    return 0


def parse_args(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("build-graph", help="build an ecspr::atom_graph directory")
    q.add_argument("--source", choices=("pairs", "gem", "gpr"), required=True)
    q.add_argument("--atom-pairs", required=True)
    q.add_argument("--direction", default=None)
    q.add_argument("--elements", nargs="*", default=None)
    q.add_argument("--weights", default=None, help="evidence_weights.parquet (source=pairs)")
    q.add_argument("--evidence-source", default="epi300")
    q.add_argument("--evidence-column", default="E_full")
    q.add_argument("--model", default=None, help="COBRA JSON or SBML (source=gem|gpr)")
    q.add_argument("--reac-xref", default=None)
    q.add_argument("--gene-set", default=None, help="active gene ids, one per line (gpr)")
    q.add_argument("--uniform-e", type=float, default=1.0)
    q.add_argument("--out-dir", required=True)
    q.set_defaults(func=cmd_build_graph)
    s = sub.add_parser("selftest")
    s.set_defaults(func=lambda a: _selftest())
    return p.parse_args(argv)


def _toy_pairs() -> pd.DataFrame:
    """Two reactions over a shared metabolite, one with two carbon transfers."""
    return pd.DataFrame([
        dict(mnxr="R1", element="C", substrate="A", product="B", sub_idx=0, prod_idx=0,
             pair_w=1.0, confidence=0.9),
        dict(mnxr="R1", element="C", substrate="A", product="B", sub_idx=1, prod_idx=1,
             pair_w=0.5, confidence=0.9),
        dict(mnxr="R2", element="C", substrate="B", product="P", sub_idx=0, prod_idx=0,
             pair_w=1.0, confidence=1.0),
        dict(mnxr="R2", element="N", substrate="B", product="Q", sub_idx=0, prod_idx=0,
             pair_w=1.0, confidence=1.0),
        dict(mnxr="R3", element="C", substrate="A", product="B", sub_idx=0, prod_idx=0,
             pair_w=1.0, confidence=1.0),          # parallel with R1's first transfer
    ])


def _selftest_pairs():
    print("[build] graph_from_pairs: weights, ratios, element filter, parallel merge")
    p = _toy_pairs()
    g = graph_from_pairs(p, "C", {"R1": 2.0, "R2": 1.0, "R3": 3.0}, {"R2": 0.01})
    e = {(g.nodes[a], g.nodes[b]): (float(g.gp[i]), float(g.gm[i]))
         for i, (a, b) in enumerate(g.edges)}
    for k in sorted(e, key=str):
        print(f"  {k}: gp={e[k][0]:.4f} gm={e[k][1]:.4f}")
    # R1 (E=2, pair_w=1) parallel with R3 (E=3, pair_w=1) on A0->B0 -> gp 5.0
    assert abs(e[(("A", 0), ("B", 0))][0] - 5.0) < 1e-12, "parallel edges must sum"
    assert abs(e[(("A", 1), ("B", 1))][0] - 1.0) < 1e-12, "pair_w must scale the edge"
    assert abs(e[(("B", 0), ("P", 0))][1] - 0.01) < 1e-12, "ratio must set gm"
    assert (("B", 0), ("Q", 0)) not in e, "the N transfer must not enter the C graph"
    assert g.meta["n_aam_gap"] == 0 and g.meta["n_reactions_used"] == 3
    # a reaction with no pair rows is the AAM gap, counted not silent
    g2 = graph_from_pairs(p, "C", {"R1": 1.0, "RX": 1.0})
    assert g2.meta["n_aam_gap"] == 1, g2.meta
    print(f"  AAM gap reported: {g2.meta['n_aam_gap']} of "
          f"{g2.meta['n_reactions_requested']} requested\n  PASS\n")
    return 0


def _selftest_schema():
    print("[build] both atom-pair schemas (int rows and comma-joined lists) agree")
    ints = pd.DataFrame([
        dict(mnxr="R1", element="C", substrate="A", product="B", sub_idx=0, prod_idx=2,
             pair_w=1.0),
        dict(mnxr="R1", element="C", substrate="A", product="B", sub_idx=1, prod_idx=3,
             pair_w=0.5),
    ])
    strs = pd.DataFrame([
        dict(mnxr="R1", element="C", substrate="A", product="B", sub_idx="0,1",
             prod_idx="2,3", pair_w="1.0,0.5"),
    ])
    a = graph_from_pairs(ints, "C", {"R1": 2.0})
    b = graph_from_pairs(strs, "C", {"R1": 2.0})
    ea = {(a.nodes[u], a.nodes[v]): float(a.gp[i]) for i, (u, v) in enumerate(a.edges)}
    eb = {(b.nodes[u], b.nodes[v]): float(b.gp[i]) for i, (u, v) in enumerate(b.edges)}
    print(f"  int schema: {ea}")
    print(f"  str schema: {eb}")
    assert ea == eb, "the two schemas must build the same graph"
    print("  PASS\n")
    return 0


def _selftest_roundtrip():
    print("[build] atom_graph directory round-trips through npz + JSON manifest")
    import tempfile
    p = _toy_pairs()
    g = graph_from_pairs(p, "C", {"R1": 2.0, "R2": 1.0, "R3": 3.0}, {"R2": 0.01})
    with tempfile.TemporaryDirectory() as d:
        write_graph_dir(d, {"C": g}, dict(source="pairs"))
        assert _elements_in(d) == ["C"], _elements_in(d)
        h = read_graph_dir(d, "C")
    print(f"  wrote/read {h.n} nodes, {h.m} edges; meta element={h.meta['element']}")
    assert h.nodes == g.nodes and h.edges == g.edges
    assert np.allclose(h.gp, g.gp) and np.allclose(h.gm, g.gm)
    assert h.meta["n_aam_gap"] == g.meta["n_aam_gap"]
    print("  PASS\n")
    return 0


def _selftest_solve_end_to_end():
    """A built graph must measure: the probe reads a share off it, and dropping a reaction
    (the caller-side LOF) moves that share."""
    print("[build] end to end: build -> solve -> per-precursor share -> caller-side LOF")
    p = pd.DataFrame([
        dict(mnxr="R1", element="C", substrate="S", product="M", sub_idx=0, prod_idx=0,
             pair_w=1.0),
        dict(mnxr="R2", element="C", substrate="M", product="P1", sub_idx=0, prod_idx=0,
             pair_w=1.0),
        dict(mnxr="R3", element="C", substrate="M", product="P2", sub_idx=0, prod_idx=0,
             pair_w=1.0),
    ])
    w = {"R1": 1.0, "R2": 2.0, "R3": 1.0}
    g = graph_from_pairs(p, "C", w)
    sol = solve(g, Terminal.metabolite(g, "S"), Terminal.merge(g, ["P1", "P2"]))
    ko = {k: v for k, v in w.items() if k != "R2"}
    gk = graph_from_pairs(p, "C", ko)
    sk = solve(gk, Terminal.metabolite(gk, "S"), Terminal.merge(gk, ["P1", "P2"]))
    print(f"  base: total={sol.total:.6f} P1={sol.share('P1'):.6f} P2={sol.share('P2'):.6f}")
    print(f"  -R2:  total={sk.total:.6f} P1={sk.share('P1'):.6f} P2={sk.share('P2'):.6f}")
    print(f"  d(total)={sk.total - sol.total:+.6f}  d(share P2)="
          f"{sk.share('P2') - sol.share('P2'):+.6f}")
    assert abs(sol.share("P1") - 2.0 / 3.0) < 1e-9, "shares must follow the conductances"
    assert sk.share("P2") == 1.0 and sk.share("P1") == 0.0
    assert sk.total < sol.total, "removing a route must lower the total"
    print("  PASS\n")
    return 0


def _selftest_gpr_polarity():
    """The wild type must be the WHOLE model, and a knockout must only ever remove.

    cobra's `GPR.eval` takes KNOCKOUTS, not active genes. Handing it the active set
    inverts the question, and the symptom is not a crash: a first run reported 532 of
    2742 reactions live for the wild type and MORE live after a knockout. Both directions
    are pinned here, on a synthetic model, so the polarity cannot silently flip back.
    """
    print("[build] GPR polarity: wild type == whole model; a knockout only removes")
    try:
        import cobra
    except ImportError:
        print("  SKIP (cobra not installed)\n")
        return 0
    m = cobra.Model("toy")
    a, b, c = (cobra.Metabolite(x) for x in ("A", "B", "C"))
    r1 = cobra.Reaction("R1"); r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("R2"); r2.add_metabolites({b: -1, c: 1})
    r3 = cobra.Reaction("R3"); r3.add_metabolites({a: -1, c: 1})
    m.add_reactions([r1, r2, r3])
    r1.gene_reaction_rule = "g1 or g2"      # isozymes: one knockout is not enough
    r2.gene_reaction_rule = "g3 and g4"     # a complex: either knockout kills it
    r3.gene_reaction_rule = ""              # ruleless: no gene to knock out
    genes = sorted(g.id for g in m.genes)
    print(f"  model genes: {genes}")

    def live(active):
        ko = {g.id for g in m.genes} - set(active)
        out = set()
        for r in m.reactions:
            rule = (r.gene_reaction_rule or "").strip()
            if not rule or bool(r.gpr.eval(ko)):
                out.add(r.id)
        return out

    wt = live(genes)
    print(f"  wild type live: {sorted(wt)}")
    assert wt == {"R1", "R2", "R3"}, f"the wild type must be the whole model: {wt}"
    one_iso = live([g for g in genes if g != "g1"])
    print(f"  -g1 (one isozyme) live: {sorted(one_iso)}")
    assert one_iso == {"R1", "R2", "R3"}, "an isozyme pair survives one knockout"
    half_complex = live([g for g in genes if g != "g3"])
    print(f"  -g3 (half a complex) live: {sorted(half_complex)}")
    assert half_complex == {"R1", "R3"}, "a complex dies when either half goes"
    nothing = live([])
    print(f"  all genes out: {sorted(nothing)}")
    assert nothing == {"R3"}, "only the ruleless reaction survives a total knockout"
    for s in (one_iso, half_complex, nothing):
        assert s <= wt, "a knockout may only ever REMOVE reactions"
    print("  PASS\n")
    return 0


def _selftest() -> int:
    print("=" * 70)
    print("ecspr_build self-tests")
    print("=" * 70)
    _selftest_pairs()
    _selftest_gpr_polarity()
    _selftest_schema()
    _selftest_roundtrip()
    _selftest_solve_end_to_end()
    print("ALL PASS")
    return 0


def main(argv=None):
    a = parse_args(argv)
    return a.func(a) or 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["selftest"]))
