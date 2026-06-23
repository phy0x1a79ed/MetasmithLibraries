"""Per-reaction evidence weight E_r from per-ORF reaction nominations.

Faithful port of scadc 04_reaction_network/10_evidence_weights.py, restricted
to the (single) fosmid source. Each ORF carries a total belief of 1.0 split
equally across the lanes that annotated it, then within a lane across
nominations by raw_score, then evenly across each nomination's MNXR fanout.
Emits E_full (all lanes) and E_dlec (dl_ec lane only) per (source, mnxr).
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
env       = model.AddRequirement(lib.GetType("containers::metabolic.env.yml"))
evidence  = model.AddRequirement(lib.GetType("metabolic::reaction_evidence"))
weights   = model.AddProduct(lib.GetType("metabolic::reaction_evidence_weights"))

_SCRIPT = r'''
import sys
import numpy as np
import pandas as pd

src_path, dst_path = sys.argv[1], sys.argv[2]


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
    nonzero = nom["sum_s"] > 0
    nom["w_n"] = np.where(nonzero, nom["s_n"] / nom["sum_s"], 1.0 / nom["n_nom"])
    nom["contrib"] = nom["w_n"] / nom["F_n"]
    out = d.merge(
        nom[["orf", "channel", "intermediate_id", "contrib"]],
        on=["orf", "channel", "intermediate_id"], how="left",
    )
    n_lanes = out.groupby("orf")["channel"].transform("nunique")
    out["contrib"] = out["contrib"] / n_lanes
    return out


def _assert_conservation(rows, label):
    per_orf = rows.groupby("orf")["contrib"].sum()
    if len(per_orf):
        worst = float((per_orf - 1.0).abs().max())
        assert worst < 1e-9, f"{label}: per-orf conservation off by {worst:.2e}"


def compute_E(df_src, label):
    rows = _nomination_contributions(df_src)
    _assert_conservation(rows, label)
    return rows.groupby("mnxr")["contrib"].sum()


ev = pd.read_parquet(src_path)
parts = []
for source, g in ev.groupby("source"):
    e_full = compute_E(g, f"{source}/full")
    dlec = g[g["channel"] == "dl_ec"]
    e_dlec = compute_E(dlec, f"{source}/dlec") if len(dlec) else pd.Series(dtype=float)
    n_orf = g.groupby("mnxr")["orf"].nunique()
    part = (
        pd.DataFrame({"E_full": e_full})
        .join(e_dlec.rename("E_dlec"), how="outer")
        .join(n_orf.rename("n_orf"), how="left")
    )
    part["E_full"] = part["E_full"].fillna(0.0)
    part["E_dlec"] = part["E_dlec"].fillna(0.0)
    part["n_orf"] = part["n_orf"].fillna(0).astype(int)
    part.insert(0, "source", source)
    part = part.reset_index().rename(columns={"index": "mnxr"})
    parts.append(part)

out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
    columns=["source", "mnxr", "E_full", "E_dlec", "n_orf"]
)
out = out[["source", "mnxr", "E_full", "E_dlec", "n_orf"]]
out.to_parquet(dst_path, index=False)
print(f"evidence_weights: {len(out)} (source,mnxr) rows")
'''


def protocol(context: ExecutionContext):
    iev = context.Input(evidence)
    iout = context.Output(weights)

    script = "evidence_weights.py"
    with open(script, "w") as f:
        f.write(_SCRIPT)
    context.ExecWithContainer(
        image=env,
        cmd=f"python {script} {iev.container} {iout.container}",
    )
    return ExecutionResult(
        manifest=[{weights: iout.local}],
        success=iout.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=evidence,
    resources=Resources(cpus=2, memory=Size.GB(8), duration=Duration(hours=1)),
)
