"""ORF homology -> per-ORF reaction nominations (reaction_evidence).

The annotation lane that connects FabFos ORFs to the metabolic graph chain:
DIAMOND blastp of the fosmid ORFs against a reference reaction DB, then project
each top subject hit to a MetaNetX reaction (MNXR) via the bundled
subject->MNXR bridge. Produces the per-ORF reaction_evidence table the
evidence_weights / build_fosmid_graph transforms consume.

This mirrors the scadc DIAMOND+Rhea lane in shape; it requires the reference
reaction DB (reactions.dmnd + bridge.tsv) to be provisioned. raw_score is the
bitscore (the evidence weighting renormalizes per ORF, so absolute scale is not
load-bearing).
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
img_dmnd  = model.AddRequirement(lib.GetType("containers::diamond.oci"))
env       = model.AddRequirement(lib.GetType("containers::metabolic.env.yml"))
orfs      = model.AddRequirement(lib.GetType("sequences::open_reading_frames"))
ref       = model.AddRequirement(lib.GetType("metabolic::reaction_reference_db"))
evidence  = model.AddProduct(lib.GetType("metabolic::reaction_evidence"))

_PROJECT = r'''
import sys
import pandas as pd

blast6, bridge_tsv, dst = sys.argv[1], sys.argv[2], sys.argv[3]

cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
hits = pd.read_csv(blast6, sep="\t", names=cols)
# top hit per ORF by bitscore
hits = hits.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")

bridge = pd.read_csv(bridge_tsv, sep="\t", names=["sseqid", "mnxr"])
merged = hits.merge(bridge, on="sseqid", how="inner")

ev = pd.DataFrame({
    "source": "fosmid",
    "orf": merged["qseqid"],
    "channel": "diamond_rhea",
    "intermediate_id": merged["sseqid"],
    "mnxr": merged["mnxr"],
    "raw_score": merged["bitscore"].astype(float),
})
ev = ev.drop_duplicates(["orf", "channel", "intermediate_id", "mnxr"])
ev.to_parquet(dst, index=False)
print(f"reaction_evidence: {len(ev)} nominations over {ev['orf'].nunique()} ORFs")
'''


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    iref = context.Input(ref)
    iev = context.Output(evidence)

    threads = context.params.get("cpus")
    tflag = "" if threads is None else f"-p {threads}"

    context.ExecWithContainer(
        image=img_dmnd,
        cmd=f"""\
            diamond blastp {tflag} \
                --query {iorfs.container} \
                --db {iref.container}/reactions.dmnd \
                --outfmt 6 \
                --max-target-seqs 5 --evalue 1e-5 \
                --out hits.b6
        """,
    )

    script = "project_reactions.py"
    with open(script, "w") as f:
        f.write(_PROJECT)
    context.ExecWithContainer(
        image=env,
        cmd=f"python {script} hits.b6 {iref.container}/bridge.tsv {iev.container}",
    )
    return ExecutionResult(
        manifest=[{evidence: iev.local}],
        success=iev.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(cpus=8, memory=Size.GB(16), duration=Duration(hours=4)),
)
