"""ECSPr belief-conservation evidence weights: evidence_table -> evidence_weights.

Thin wrapper around `resources/lib/ecspr_network.py weights` (verified byte-exact
vs scadc evidence_weights.parquet). Allocates each ORF's belief 1.0 across its
lanes -> raw_score -> 1/fanout over MNXR, giving per-(source,orf,mnxr) E_full.
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
evidence = model.AddRequirement(lib.GetType("functional_annotation::evidence_table"), parents={exp})
env      = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
net      = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver   = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
weights  = model.AddProduct(lib.GetType("ecspr::evidence_weights"))

def protocol(context: ExecutionContext):
    iev  = context.Input(evidence)
    inet = context.Input(net)
    owet = context.Output(weights)
    context.ExecWithContainer(
        image=env,
        cmd=f"python {inet.container} weights --evidence {iev.container} --out {owet.container}",
    )
    return ExecutionResult(
        manifest=[{weights: owet.local}],
        success=owet.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(minutes=30)),
)
