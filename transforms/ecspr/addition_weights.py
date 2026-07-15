"""ECSPr fosmid addition maps: evidence_table -> addition_weights pickle.

Thin wrapper around `resources/lib/ecspr_network.py addition-weights`
(--source fosmid --unit contig): the per-fosmid {contig: {mnxr: E}} reaction
maps injected onto the host base graph in the SMW solve.
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
evidence = model.AddRequirement(lib.GetType("functional_annotation::evidence_table"), parents={exp})
env      = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
net      = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver   = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
addition = model.AddProduct(lib.GetType("ecspr::addition_weights"))

def protocol(context: ExecutionContext):
    iev  = context.Input(evidence)
    inet = context.Input(net)
    oadd = context.Output(addition)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {inet.container} addition-weights \
            --evidence {iev.container} --source fosmid --unit contig \
            --out {oadd.container}""",
    )
    return ExecutionResult(
        manifest=[{addition: oadd.local}],
        success=oadd.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(minutes=30)),
)
