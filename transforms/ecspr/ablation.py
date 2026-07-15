"""ECSPr per-gene ablation: evidence + addition + base_graphs + bipartite + axes
-> ablation_importance (LOO + gene-solo, BOTH metrics).

Thin wrapper around `resources/lib/ecspr_ablation.py` (verified vs scadc
reff_ablation_importance.tsv, importance ~4.6e-9). For each fosmid on each testable
axis and each ORF: full-insert delta, leave-one-out ("minus each gene"), and
gene-solo ("each gene alone") -- for BOTH effective resistance (delta_reff) and
effective conductance (delta_ieff). Reads axes_testable.json from the base dir.

device/dtype arrive as a STAGED `ecspr::compute_profile` input, never via
context.params -- see transforms/ecspr/solve.py for why that distinction is
load-bearing rather than stylistic.
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
profile   = model.AddRequirement(lib.GetType("ecspr::compute_profile"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
abl       = model.AddRequirement(lib.GetType("lib::ecspr_ablation.py"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
importance = model.AddProduct(lib.GetType("ecspr::ablation_importance"))

def read_compute_profile(path):
    """`key: value` per line -> (device, dtype). Fails loudly on an unknown value:
    a typo'd device must not silently degrade to cpu after hours of queueing.

    Deliberately duplicated from transforms/ecspr/solve.py rather than imported.
    A transform's own directory IS on sys.path at load time, so `from solve import
    ...` resolves -- but it also EXECUTES solve.py's module body, registering solve's
    TransformInstance as a side effect of loading this one. Ten duplicated lines beat
    a planner that mis-registers transforms in a way that reads as a design problem.
    """
    prof = {}
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        k, _, v = line.partition(":")
        prof[k.strip()] = v.strip()
    device, dtype = prof.get("device", "cpu"), prof.get("dtype", "float64")
    assert device in ("cpu", "cuda"), f"compute_profile: bad device {device!r}"
    assert dtype in ("float64", "float32"), f"compute_profile: bad dtype {dtype!r}"
    return device, dtype

def protocol(context: ExecutionContext):
    iev   = context.Input(evidence)
    iadd  = context.Input(addition)
    ibase = context.Input(bases)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    iabl  = context.Input(abl)
    iprof = context.Input(profile)
    oimp  = context.Output(importance)
    device, dtype = read_compute_profile(iprof.local)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {iabl.container} \
            --evidence {iev.container} --addition {iadd.container} \
            --base-dir {ibase.container} --bipartite-dir {ibip.container} \
            --axes {iaxes.container} --testable {ibase.container}/axes_testable.json \
            --out {oimp.container} --device {device} --dtype {dtype}""",
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
