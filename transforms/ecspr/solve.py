"""ECSPr solve: addition + base_graphs + bipartite + axes -> reff + ieff reports.

Thin wrapper around `resources/lib/ecspr_network.py solve`. One pass per fosmid
yields BOTH effective resistance (delta_reff) and effective conductance
(delta_ieff = 1/r_aug - 1/r_base) on every axis. Reads axes_testable.json from the
base_graphs dir.

ON THE CLAIM THIS DOCSTRING USED TO MAKE
----------------------------------------
It said "verified vs scadc reff_axes_report.tsv, max|abs| 1.4e-8". That sentence was
the only thing standing in for a solver gate, it had never been executed, and 1.4e-8
is an order of magnitude LOOSER than canon.PARITY_TOL -- the tolerance the same repo
calls non-negotiable. Worse, the referent was wrong: measured against a from-scratch
rebuild, this engine agrees to ~1e-13 while that incumbent table agrees only to
~1e-7. The table is the inaccurate artifact, so "verified against it" was aiming at
the wrong target.

There is now a real gate, and it is not in this docstring:

    python parity/run_solver_parity.py     # in the scadc fabfos methods path

It scores BOTH solver paths against a rebuild that shares none of their math.

WHY THE DEFAULT PATH IS `direct` AND NOT WOODBURY
-------------------------------------------------
Woodbury pays only when Z = L^-1 is a precomputed DENSE inverse, so that Z_ee is a
gather -- true on GPU, false on CPU, where it costs p sparse solves per (axis,
fosmid). Measured on carbon: 88.8s vs 1.36s for 8 fosmids x 40 axes. Over the basis,
~45 min -> ~63s. `--path woodbury` is kept reachable because it is independently
verified and therefore serves as this path's referent.

device/dtype arrive as a STAGED `ecspr::compute_profile` input, never via
context.params. params are populated only from the runtime resources line
(cpus/memory/attempt), so `context.params.get("device")` was ALWAYS its default --
the GPU path was unreachable through the driver, which is why an earlier
end-to-end run bypassed metasmith with raw sbatch. params also never enter the
task hash, so two runs intending different knobs would collide on one cache entry.
As an input the profile is content-hashed: a GPU run and a CPU run of this step
produce different cache keys, which is the property that makes the whole thing safe.

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
