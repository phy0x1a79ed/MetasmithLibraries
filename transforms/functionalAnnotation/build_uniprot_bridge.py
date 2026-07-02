"""UniProt -> MNXR bridge via Rhea DR (lane-3 projection).

Thin wrapper around `resources/lib/fabfos_evidence.py build-uniprot-bridge`. Joins
rhea2uniprot{,_trembl} against reac_xref's rhea: rows to map each UniProt accession
to its MNXR(s). Pure reference transform (MetaNetX + Rhea). The TrEMBL table is the
big one (~hundreds of MB); reviewed-only still works if it is not staged.
"""
from metasmith.python_api import *

lib        = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model      = Transform()
reac_xref  = model.AddRequirement(lib.GetType("functional_annotation::metanetx_reac_xref"))
rhea_sw    = model.AddRequirement(lib.GetType("functional_annotation::rhea2uniprot"))
rhea_tr    = model.AddRequirement(lib.GetType("functional_annotation::rhea2uniprot_trembl"))
env        = model.AddRequirement(lib.GetType("envs::ecspr.condaenv"))
ev         = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
bridge     = model.AddProduct(lib.GetType("functional_annotation::uniprot_to_mnxr"))

def protocol(context: ExecutionContext):
    ixref = context.Input(reac_xref)
    isw   = context.Input(rhea_sw)
    itr   = context.Input(rhea_tr)
    iev   = context.Input(ev)
    obr   = context.Output(bridge)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {iev.container} build-uniprot-bridge \
            --reac-xref {ixref.container} \
            --rhea-swiss {isw.container} \
            --rhea-trembl {itr.container} \
            --out {obr.container}""",
    )
    return ExecutionResult(
        manifest=[{bridge: obr.local}],
        success=obr.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=reac_xref,
    resources=Resources(cpus=1, memory=Size.GB(16), duration=Duration(hours=2)),
)
