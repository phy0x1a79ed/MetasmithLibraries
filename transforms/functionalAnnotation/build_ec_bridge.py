"""EC(level-4) -> MNXR bridge from MetaNetX reac_prop classifs (lane-2 projection).

Thin wrapper around `resources/lib/fabfos_evidence.py build-ec-bridge`. Parses the
`;`-separated EC classifs column of reac_prop.tsv into a long ec->mnxr fan-out.
Pure reference transform (MetaNetX only) -- runs anywhere pandas is present.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
reac_prop = model.AddRequirement(lib.GetType("ecspr::metanetx_reac_prop"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
ev        = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
bridge    = model.AddProduct(lib.GetType("functional_annotation::ec_to_mnxr"))

def protocol(context: ExecutionContext):
    irp  = context.Input(reac_prop)
    iev  = context.Input(ev)
    obr  = context.Output(bridge)
    context.ExecWithContainer(
        image=env,
        cmd=f"python {iev.container} build-ec-bridge --reac-prop {irp.container} --out {obr.container}",
    )
    return ExecutionResult(
        manifest=[{bridge: obr.local}],
        success=obr.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=reac_prop,
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(minutes=20)),
)
