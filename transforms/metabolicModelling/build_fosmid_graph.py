"""Per-fosmid element-weighted bipartite metabolic graphs.

The endpoint of the FabFos -> ECSPr prerequisite chain: for each fosmid and
each element (C/N/S/P), overlay the fosmid's evidence-weighted reactions onto
the host base graph, producing exactly the graphs the ECSPr axis/conductance
step consumes. This transform STOPS here — no current-flow betweenness, delta
conductance, null model, or significance is computed.

Ports scadc 04_reaction_network/13_base_graphs.py (per-contig addition weights)
and _smw.build_ar2m (the evidence-weighted overlay: conductance w_X * E_r, kept
to metabolites already in base, in-base reactions entering as parallel
`rxn_reinf` copies). Base graphs and the universe element bipartite are supplied
as bundled reference inputs rather than recomputed.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
env       = model.AddRequirement(lib.GetType("containers::metabolic.env.yml"))
evidence  = model.AddRequirement(lib.GetType("metabolic::reaction_evidence"))
universe  = model.AddRequirement(lib.GetType("metabolic::element_bipartite"))
bases     = model.AddRequirement(lib.GetType("metabolic::base_graphs"))
add_w     = model.AddProduct(lib.GetType("metabolic::fosmid_addition_weights"))
graphs    = model.AddProduct(lib.GetType("metabolic::fosmid_bipartite_graph"))

_SCRIPT = r'''
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

ev_path     = Path(sys.argv[1])     # reaction_evidence parquet (fosmid source)
universe_dir = Path(sys.argv[2])    # mnx_bipartite_{C,N,S,P}.pkl
base_dir    = Path(sys.argv[3])     # base_{C,N,S,P}.pkl
add_out     = Path(sys.argv[4])     # fosmid_addition_weights.pkl
graph_dir   = Path(sys.argv[5])     # output dir of fosmid_graph_{X}.pkl
graph_dir.mkdir(parents=True, exist_ok=True)

ELEMENTS = ["C", "N", "S", "P"]


def _nomination_contributions(df):
    d = df.drop_duplicates(["orf", "channel", "intermediate_id", "mnxr"]).copy()
    nom = (
        d.groupby(["orf", "channel", "intermediate_id"], sort=False)
        .agg(s_n=("raw_score", "max"), F_n=("mnxr", "nunique"))
        .reset_index()
    )
    grp = nom.groupby(["orf", "channel"], sort=False)
    nom["sum_s"] = grp["s_n"].transform("sum")
    nom["n_nom"] = grp["intermediate_id"].transform("count")
    nz = nom["sum_s"] > 0
    nom["w_n"] = np.where(nz, nom["s_n"] / nom["sum_s"], 1.0 / nom["n_nom"])
    nom["contrib"] = nom["w_n"] / nom["F_n"]
    out = d.merge(
        nom[["orf", "channel", "intermediate_id", "contrib"]],
        on=["orf", "channel", "intermediate_id"], how="left",
    )
    n_lanes = out.groupby("orf")["channel"].transform("nunique")
    out["contrib"] = out["contrib"] / n_lanes
    return out


def per_contig_weights(df):
    """{contig: {mnxr: E}} — per-fosmid additive evidence weights."""
    rows = _nomination_contributions(df).copy()
    rows["contig"] = rows["orf"].str.rsplit("_", n=1).str[0]
    agg = rows.groupby(["contig", "mnxr"], sort=False)["contrib"].sum()
    out = defaultdict(dict)
    for (contig, mnxr), c in agg.items():
        if c > 0:
            out[contig][mnxr] = float(c)
    return dict(out)


def build_ar2m(weights, G_X, base_nodes, wkey):
    """Evidence-weighted addition map; ports _smw.build_ar2m (reinforce=True)."""
    ar2m = {}
    for mnxr, er in weights.items():
        if er <= 0:
            continue
        canon = ("rxn", mnxr)
        if canon not in G_X:
            continue
        node_id = ("rxn_reinf", mnxr) if canon in base_nodes else canon
        entries = []
        for _, m, d in G_X.edges(canon, data=True):
            if m[0] != "met" or m not in base_nodes:
                continue
            w = d.get(wkey, 0.0)
            if w > 0:
                entries.append((m, float(w) * float(er)))
        if entries:
            ar2m[node_id] = entries
    return ar2m


ev = pd.read_parquet(ev_path)
fos = ev[ev["source"] == "fosmid"] if "source" in ev.columns else ev
contig_w = per_contig_weights(fos)
with open(add_out, "wb") as fh:
    pickle.dump(contig_w, fh, protocol=pickle.HIGHEST_PROTOCOL)
print(f"fosmid_addition_weights: {len(contig_w)} contigs")

for X in ELEMENTS:
    bpkl = base_dir / f"base_{X}.pkl"
    upkl = universe_dir / f"mnx_bipartite_{X}.pkl"
    if not (bpkl.exists() and upkl.exists()):
        print(f"  [{X}] missing base/universe pickle; skipping")
        continue
    with open(bpkl, "rb") as fh:
        base = pickle.load(fh)
    with open(upkl, "rb") as fh:
        G_X = pickle.load(fh)
    base_nodes = set(base.nodes)
    wkey = f"w_{X}"
    per_fosmid = {}
    for contig, weights in contig_w.items():
        ar2m = build_ar2m(weights, G_X, base_nodes, wkey)
        G = base.copy()
        for node_id, entries in ar2m.items():
            for m, c in entries:
                G.add_edge(node_id, m, **{wkey: c})
        per_fosmid[contig] = G
    with open(graph_dir / f"fosmid_graph_{X}.pkl", "wb") as fh:
        pickle.dump(per_fosmid, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  [{X}] {len(per_fosmid)} per-fosmid graphs "
          f"(base {base.number_of_nodes()} nodes)")
'''


def protocol(context: ExecutionContext):
    iev = context.Input(evidence)
    iuni = context.Input(universe)
    ibase = context.Input(bases)
    iaddw = context.Output(add_w)
    igraphs = context.Output(graphs)

    script = "build_fosmid_graph.py"
    with open(script, "w") as f:
        f.write(_SCRIPT)
    context.ExecWithContainer(
        image=env,
        cmd=(
            f"python {script} {iev.container} {iuni.container} {ibase.container} "
            f"{iaddw.container} {igraphs.container}"
        ),
    )
    return ExecutionResult(
        manifest=[{add_w: iaddw.local, graphs: igraphs.local}],
        success=iaddw.local.exists() and igraphs.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=evidence,
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=2)),
)
