"""ECSPr DIRECTED solve: addition + base_graphs + bipartite + axes + roles + direction
-> reff + ieff reports (the rectified-network / diode analog of the undirected solve).

Thin wrapper around `resources/lib/ecspr_network.py solve-directed`. Orients the
reference base per reaction ROLES (reac_prop.tsv -> MNXR->(S,P)) and the direction
ensemble's g_rev/g_fwd RATIO, then solves every (fosmid, axis) as a rectified-network
C_eff through the engine's directed primitive (ecspr_directed). Emits the SAME
reff/ieff schema as the undirected solve, so the significance scorer consumes it
UNCHANGED -- the directed delta_ieff flows through the canonical lane, no new column,
no new lane.

WHY THIS IS A SEPARATE DOMAIN (ecsprDirected), NOT A SECOND transforms/ecspr/*.py
--------------------------------------------------------------------------------
Both this and the undirected `ecspr/solve.py` produce `ecspr::reff_axes_report` +
`ecspr::ieff_axes_report`. A loaded domain is selectable structurally by the planner,
so if both lived in one loaded domain the planner would face two producers of the same
type and pick one by tiebreak -- exactly the silent-lane failure canon exists to
prevent. Isolating the directed solve in its own domain makes producer selection a
DELIBERATE choice of which domain the experiment loads (the same precedent ecsprNetA
set for the curated-GEM base builder), never a planner coin-flip. See
main/fabfos/experiments and MIGRATION.md 2026-07-18.

WHY THERE IS NO compute_profile HERE
------------------------------------
The directed solve is a CPU semismooth-Newton solve: `solve-directed` has NO
--device/--dtype argument (unlike the undirected `solve`), so there is nothing for a
compute_profile to carry -- staging one would be a hashed input the CLI ignores. The
only numeric knob is Newton --tol, left at the engine default (committed there, not
here).

WHY ecspr_directed.py IS REQUIRED THOUGH THE PROTOCOL NEVER READS IT
-------------------------------------------------------------------
`ecspr_network.py` imports `ecspr_directed` at module top, so solve-directed cannot
import without it staged beside ecspr_network.py in the container -- exactly as the
undirected solve requires ecspr_solver.py without context.Input-ing it. Requiring the
lib resource stages it into the shared container dir; the sibling import then resolves.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
addition  = model.AddRequirement(lib.GetType("ecspr::addition_weights"), parents={exp})
bases     = model.AddRequirement(lib.GetType("ecspr::base_graphs"), parents={exp})
bipartite = model.AddRequirement(lib.GetType("ecspr::mnx_bipartite"))
axes      = model.AddRequirement(lib.GetType("ecspr::biomass_axes"))
roles     = model.AddRequirement(lib.GetType("ecspr::reaction_roles"))
direction = model.AddRequirement(lib.GetType("ecspr::direction_ratios"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
directed  = model.AddRequirement(lib.GetType("lib::ecspr_directed.py"))
reff      = model.AddProduct(lib.GetType("ecspr::reff_axes_report"))
ieff      = model.AddProduct(lib.GetType("ecspr::ieff_axes_report"))

def protocol(context: ExecutionContext):
    iadd  = context.Input(addition)
    ibase = context.Input(bases)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    iroles = context.Input(roles)
    idir  = context.Input(direction)
    inet  = context.Input(net)
    oreff = context.Output(reff)
    oieff = context.Output(ieff)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {inet.container} solve-directed \
            --addition {iadd.container} --base-dir {ibase.container} \
            --bipartite-dir {ibip.container} --axes {iaxes.container} \
            --testable {ibase.container}/axes_testable.json \
            --roles {iroles.container} --direction {idir.container} \
            --out-reff {oreff.container} --out-ieff {oieff.container}""",
    )
    return ExecutionResult(
        manifest=[{reff: oreff.local, ieff: oieff.local}],
        success=oreff.local.exists() and oieff.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    # cpus=1, not 4. `solve-directed` takes no thread flag and no --device, and
    # the solve factorizes through scipy's SuperLU, which is serial -- so the
    # extra three cores were reserved and then idled. Worse, the engine pins
    # OMP/OPENBLAS/MKL/NUMEXPR to 1 thread on purpose (oversubscribing measured
    # >10x slower), so even the BLAS underneath will not use them. Parallelism
    # here comes from running many solves at once, not from widening one.
    resources=Resources(cpus=1, memory=Size.GB(32), duration=Duration(hours=12)),
)
