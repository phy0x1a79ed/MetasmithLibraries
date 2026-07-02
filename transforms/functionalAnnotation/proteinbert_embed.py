"""Lane 4 (part a) -- per-ORF ProteinBERT embeddings for the fosmid ORFs.

Fresh ProteinBERT inference (mamba, fabfos-ml env, GPU) over the reference-insert
ORFs -> protein_embeddings (.npy, row i = ORF i) with a sibling .index.csv. These
are the query embeddings the lane-4 embed-transfer kNN projects labels onto.
ProteinBERT is the primary backbone (raw embeddings won its dark-regime validation);
ESM-C 600M is the optional comparison backbone (esmc_model / esmc.condaenv).

GATED: needs the ProteinBERT weights (proteinbert_model dir) + fabfos-ml env
(torch) and a GPU to be practical. The `proteinbert` CLI must be on PATH.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
exp   = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
orfs  = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
mdl   = model.AddRequirement(lib.GetType("functional_annotation::proteinbert_model"))
env   = model.AddRequirement(lib.GetType("envs::ml.condaenv"))
emb   = model.AddProduct(lib.GetType("functional_annotation::protein_embeddings"))

def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    imdl  = context.Input(mdl)
    oemb  = context.Output(emb)
    device = context.params.get("device", "cuda")
    context.ExecWithContainer(
        image=env,
        cmd=f"""proteinbert embed \
            --model {imdl.container} \
            --fasta {iorfs.container} \
            --device {device} \
            --out {oemb.container}""",
    )
    return ExecutionResult(
        manifest=[{emb: oemb.local}],
        success=oemb.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=4)),
)
