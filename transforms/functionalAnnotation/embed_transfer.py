"""Lane 4 (part b) -- embedding-transfer kNN EC/reaction labels -> embed_transfer.

Thin wrapper around `resources/lib/fabfos_embed_transfer.py apply`. Distance-weighted
kNN vote over L2-normalized ProteinBERT embeddings transfers MNXR labels from the
labelled reference pool to the dark fosmid ORFs (those with zero evidence from the
other 3 lanes). ProteinBERT raw = pbert_transfer (floor 0.20); projected ESM-C =
esmc_transfer (floor 0.10). projection_via='embedding_knn'. Exploratory lane; fold
into the evidence table via the compiler's --embed once validated.

GATED: needs the staged reference_label_pool (labelled reference embeddings +
orf_index + trained projector, built by fabfos_embed_transfer.py build-reference-pool
+ train-projector over the fresh fosmid protein_embeddings + the reused reference
embeddings), plus fabfos-ml env + GPU. The reference embedding pool is a reused
REFERENCE artifact (not the metaG null).
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
exp   = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
emb   = model.AddRequirement(lib.GetType("functional_annotation::protein_embeddings"), parents={exp})
pool  = model.AddRequirement(lib.GetType("functional_annotation::reference_label_pool"))
env   = model.AddRequirement(lib.GetType("envs::ml.condaenv"))
et    = model.AddRequirement(lib.GetType("lib::fabfos_embed_transfer.py"))
out   = model.AddProduct(lib.GetType("functional_annotation::embed_transfer"))

def protocol(context: ExecutionContext):
    ipool = context.Input(pool)
    iet   = context.Input(et)
    oout  = context.Output(out)
    device = context.params.get("device", "cuda")
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {iet.container} apply \
            --in-dir {ipool.container} \
            --query-source fosmid \
            --device {device} \
            --out {oout.container}""",
    )
    return ExecutionResult(
        manifest=[{out: oout.local}],
        success=oout.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=4)),
)
