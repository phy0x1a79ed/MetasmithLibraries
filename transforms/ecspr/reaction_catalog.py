"""ECSPr reaction catalog: evidence_table + evidence_weights + MetaNetX -> 5 parquets.

Thin wrapper around `resources/lib/ecspr_catalog.py build` (selftest-verified).
Synthesizes reactions / metabolites / rxn_edges / orf_reactions / reaction_redundancy
into one directory. reac_prop + chem_prop are reused MetaNetX reference tables.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
evidence  = model.AddRequirement(lib.GetType("functional_annotation::evidence_table"), parents={exp})
weights   = model.AddRequirement(lib.GetType("ecspr::evidence_weights"), parents={exp})
reac_prop = model.AddRequirement(lib.GetType("ecspr::metanetx_reac_prop"))
chem_prop = model.AddRequirement(lib.GetType("ecspr::metanetx_chem_prop"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
cat       = model.AddRequirement(lib.GetType("lib::ecspr_catalog.py"))
net       = model.AddRequirement(lib.GetType("lib::ecspr_network.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
catalog   = model.AddProduct(lib.GetType("ecspr::reaction_catalog"))

def protocol(context: ExecutionContext):
    iev  = context.Input(evidence)
    iwet = context.Input(weights)
    irp  = context.Input(reac_prop)
    icp  = context.Input(chem_prop)
    icat = context.Input(cat)
    ocat = context.Output(catalog)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {icat.container} build \
            --evidence {iev.container} --weights {iwet.container} \
            --reac-prop {irp.container} --chem-prop {icp.container} \
            --out {ocat.container}""",
    )
    return ExecutionResult(
        manifest=[{catalog: ocat.local}],
        success=ocat.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=1)),
)
