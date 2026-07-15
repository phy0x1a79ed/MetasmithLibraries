"""ECSPr SMW solve: addition + base_graphs + bipartite + axes -> reff + ieff reports.

Thin wrapper around `resources/lib/ecspr_network.py solve` (verified vs scadc
reff_axes_report.tsv, max|abs| 1.4e-8). One SMW/Woodbury pass per (fosmid, axis)
yields BOTH effective resistance (delta_reff) and effective conductance
(delta_ieff = 1/r_aug - 1/r_base). Reads axes_testable.json from the base_graphs
dir.

device/dtype arrive as a STAGED `ecspr::compute_profile` input, never via
context.params. params are populated only from the runtime resources line
(cpus/memory/attempt), so `context.params.get("device")` was ALWAYS its default --
the GPU path was unreachable through the driver, which is why an earlier
end-to-end run bypassed metasmith with raw sbatch. params also never enter the
task hash, so two runs intending different knobs would collide on one cache entry.
As an input the profile is content-hashed: a GPU run and a CPU run of this step
produce different cache keys, which is the property that makes the whole thing safe.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
addition  = model.AddRequirement(lib.GetType("ecspr::addition_weights"), parents={exp})
bases     = model.AddRequirement(lib.GetType("ecspr::base_graphs"), parents={exp})
bipartite = model.AddRequirement(lib.GetType("ecspr::mnx_bipartite"))
axes      = model.AddRequirement(lib.GetType("ecspr::biomass_axes"))
profile   = model.AddRequirement(lib.GetType("ecspr::compute_profile"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
reff      = model.AddProduct(lib.GetType("ecspr::reff_axes_report"))
ieff      = model.AddProduct(lib.GetType("ecspr::ieff_axes_report"))

def read_compute_profile(path):
    """`key: value` per line -> (device, dtype). Fails loudly on an unknown value:
    a typo'd device must not silently degrade to cpu after hours of queueing."""
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
    iadd  = context.Input(addition)
    ibase = context.Input(bases)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    inet  = context.Input(net)
    iprof = context.Input(profile)
    oreff = context.Output(reff)
    oieff = context.Output(ieff)
    device, dtype = read_compute_profile(iprof.local)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {inet.container} solve \
            --addition {iadd.container} --base-dir {ibase.container} \
            --bipartite-dir {ibip.container} --axes {iaxes.container} \
            --testable {ibase.container}/axes_testable.json \
            --out-reff {oreff.container} --out-ieff {oieff.container} \
            --device {device} --dtype {dtype}""",
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
