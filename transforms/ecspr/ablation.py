"""ECSPr per-gene ablation: evidence + addition + base_graphs + bipartite + axes
-> ablation_importance (LOO + gene-solo, BOTH metrics).

Thin wrapper around `resources/lib/ecspr_ablation.py` (verified vs scadc
reff_ablation_importance.tsv, importance ~4.6e-9). For each fosmid on each testable
axis and each ORF: full-insert delta, leave-one-out ("minus each gene"), and
gene-solo ("each gene alone") -- for BOTH effective resistance (delta_reff) and
effective conductance (delta_ieff). Reads axes_testable.json from the base dir.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
evidence  = model.AddRequirement(lib.GetType("functional_annotation::evidence_table"), parents={exp})
addition  = model.AddRequirement(lib.GetType("ecspr::addition_weights"), parents={exp})
bases     = model.AddRequirement(lib.GetType("ecspr::base_graphs"), parents={exp})
bipartite = model.AddRequirement(lib.GetType("ecspr::mnx_bipartite"))
axes      = model.AddRequirement(lib.GetType("ecspr::biomass_axes"))
env       = model.AddRequirement(lib.GetType("envs::ecspr.condaenv"))
abl       = model.AddRequirement(lib.GetType("lib::ecspr_ablation.py"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
importance = model.AddProduct(lib.GetType("ecspr::ablation_importance"))

def protocol(context: ExecutionContext):
    iev   = context.Input(evidence)
    iadd  = context.Input(addition)
    ibase = context.Input(bases)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    iabl  = context.Input(abl)
    oimp  = context.Output(importance)
    device = context.params.get("device", "cpu")
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {iabl.container} \
            --evidence {iev.container} --addition {iadd.container} \
            --base-dir {ibase.container} --bipartite-dir {ibip.container} \
            --axes {iaxes.container} --testable {ibase.container}/axes_testable.json \
            --out {oimp.container} --device {device}""",
    )
    return ExecutionResult(
        manifest=[{importance: oimp.local}],
        success=oimp.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=12)),
)
