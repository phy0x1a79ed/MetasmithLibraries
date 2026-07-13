"""ECSPr Network A base graphs (direct curated GEM) -> base_graphs dir.

Thin wrapper around `resources/lib/ecspr_network.py gem-base-graphs`. This is the
CURATED-GEM alternative to the annotation-derived `base_graphs.py`: instead of
inducing the universe on the 4-lane fosmid evidence table, it crosswalks a curated
GEM's reactome (iML1515 for K-12, iECDH10B for DH10B/forsberg) to current MetaNetX
MNXR (via reac_xref, preferring atom-mapped universe candidates) and induces the
same atom-mapped element universe with UNIFORM conductance E=1.0 (w_X = pure
atom-transit count; no annotation over-nomination, no evidence weighting).

Reuse-only AAM: every GEM MNXR carrying a w_X>0 edge in the universe is an
atom-mapped base node -- zero new DL mapping. Produces the identical `base_graphs`
directory (base_{C,N,S,P}.pkl + axes_testable.json) so the downstream solve /
ablation / significance transforms consume Network A and Network B identically.
mnx_bipartite + reac_xref are reused reference tables (staged inputs).
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
gem       = model.AddRequirement(lib.GetType("ecspr::curated_gem"))
reac_xref = model.AddRequirement(lib.GetType("ecspr::metanetx_reac_xref"))
bipartite = model.AddRequirement(lib.GetType("ecspr::mnx_bipartite"))
axes      = model.AddRequirement(lib.GetType("ecspr::biomass_axes"))
env       = model.AddRequirement(lib.GetType("envs::ecspr.condaenv"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
bases     = model.AddProduct(lib.GetType("ecspr::base_graphs"))

def protocol(context: ExecutionContext):
    igem  = context.Input(gem)
    ixref = context.Input(reac_xref)
    ibip  = context.Input(bipartite)
    iaxes = context.Input(axes)
    inet  = context.Input(net)
    obase = context.Output(bases)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {inet.container} gem-base-graphs \
            --model {igem.container} --reac-xref {ixref.container} \
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
    group_by=gem,
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=2)),
)
