"""ECSPr SMW solve: addition + base_graphs + bipartite + axes -> reff + ieff reports.

Thin wrapper around `resources/lib/ecspr_network.py solve` (verified vs scadc
reff_axes_report.tsv, max|abs| 1.4e-8). One SMW/Woodbury pass per (fosmid, axis)
yields BOTH effective resistance (delta_reff) and effective conductance
(delta_ieff = 1/r_aug - 1/r_base). Reads axes_testable.json from the base_graphs
dir. CPU (float64); set param device=cuda for GPU float32.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
addition  = model.AddRequirement(lib.GetType("ecspr::addition_weights"), parents={exp})
bases     = model.AddRequirement(lib.GetType("ecspr::base_graphs"), parents={exp})
bipartite = model.AddRequirement(lib.GetType("ecspr::mnx_bipartite"))
axes      = model.AddRequirement(lib.GetType("ecspr::biomass_axes"))
env       = model.AddRequirement(lib.GetType("envs::ecspr.condaenv"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
reff      = model.AddProduct(lib.GetType("ecspr::reff_axes_report"))
ieff      = model.AddProduct(lib.GetType("ecspr::ieff_axes_report"))

def protocol(context: ExecutionContext):
    iadd  = context.Input(addition)
    ibase = context.Input(bases)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    inet  = context.Input(net)
    oreff = context.Output(reff)
    oieff = context.Output(ieff)
    device = context.params.get("device", "cpu")
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {inet.container} solve \
            --addition {iadd.container} --base-dir {ibase.container} \
            --bipartite-dir {ibip.container} --axes {iaxes.container} \
            --testable {ibase.container}/axes_testable.json \
            --out-reff {oreff.container} --out-ieff {oieff.container} \
            --device {device}""",
    )
    return ExecutionResult(
        manifest=[{reff: oreff.local, ieff: oieff.local}],
        success=oreff.local.exists() and oieff.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=6)),
)
