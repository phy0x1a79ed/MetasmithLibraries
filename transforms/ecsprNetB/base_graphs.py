"""ECSPr host base graphs: evidence + weights + bipartite + axes -> base_graphs dir.

Thin wrapper around `resources/lib/ecspr_network.py base-graphs`: builds the
per-element epi300 host reaction-network graph (base_{X}.pkl) and, given the
biomass axes, writes axes_testable.json (axes whose endpoints are in the element
LCC). mnx_bipartite is the reused per-element bipartite reference directory.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
evidence  = model.AddRequirement(lib.GetType("functional_annotation::evidence_table"), parents={exp})
weights   = model.AddRequirement(lib.GetType("ecspr::evidence_weights"), parents={exp})
bipartite = model.AddRequirement(lib.GetType("ecspr::mnx_bipartite"))
axes      = model.AddRequirement(lib.GetType("ecspr::biomass_axes"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
bases     = model.AddProduct(lib.GetType("ecspr::base_graphs"))

def protocol(context: ExecutionContext):
    iev   = context.Input(evidence)
    iwet  = context.Input(weights)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    inet  = context.Input(net)
    obase = context.Output(bases)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {inet.container} base-graphs \
            --evidence {iev.container} --weights {iwet.container} \
            --bipartite-dir {ibip.container} --axes {iaxes.container} \
            --out-dir {obase.container}""",
    )
    return ExecutionResult(
        manifest=[{bases: obase.local}],
        success=obase.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=2)),
)
