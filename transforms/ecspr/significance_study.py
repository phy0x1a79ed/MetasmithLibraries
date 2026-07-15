"""ECSPr null-family STUDY: reports + staged null + a staged family spec
-> {reff,ieff}_significance_study.

A new null family is ONE file: add it to FAMILIES in resources/lib/ecspr_null_study.py
and stage a spec naming it. Nothing else changes -- not this transform, not the
canonical scorer, not the DAG shape.

The family arrives as the STAGED, CONTENT-HASHED `ecspr::null_model_spec` input, so
two families produce different cache keys and cannot overwrite each other's results.
That is the whole point: the incumbent's draws file carried no K in its name, so a
K=10000 run silently overwrote the K=1000 set in place, mid-run, with the file's
mtime moving BACKWARDS. It was recovered only because the seed happened to be
deterministic. Here that class of bug is not mitigated, it is unrepresentable --
outputs are addressed by the hash of their inputs.

This produces a STUDY type. The canonical `ecspr::{lane}_significance` keeps exactly
one producer (the full-mixture SF); no staged input can change what it means. See
resources/lib/ecspr_null_study.py for why that asymmetry is deliberate.

Dev loop -- iterate a family standalone, no solver and no nextflow, but through the
same ExecuteStep the DAG uses, so what you tested is what runs:

    from metasmith.python_api import RunTransform
    RunTransform(
        transform_lib=<lib>/transforms/ecspr,
        transform="significance_study.py",
        inputs=[("ecspr::ieff_axes_report", ...), ("ecspr::null_model_spec", ...), ...],
        work_dir=<scratch>,
    )
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
reff     = model.AddRequirement(lib.GetType("ecspr::reff_axes_report"), parents={exp})
ieff     = model.AddRequirement(lib.GetType("ecspr::ieff_axes_report"), parents={exp})
proteins = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
nulls    = model.AddRequirement(lib.GetType("ecspr::frozen_null"))
famspec  = model.AddRequirement(lib.GetType("ecspr::null_model_spec"))
env      = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
study    = model.AddRequirement(lib.GetType("lib::ecspr_null_study.py"))
sig      = model.AddRequirement(lib.GetType("lib::ecspr_significance.py"))
reff_out = model.AddProduct(lib.GetType("ecspr::reff_significance_study"))
ieff_out = model.AddProduct(lib.GetType("ecspr::ieff_significance_study"))

def protocol(context: ExecutionContext):
    ireff = context.Input(reff)
    iieff = context.Input(ieff)
    ifaa  = context.Input(proteins)
    inul  = context.Input(nulls)
    ifam  = context.Input(famspec)
    istud = context.Input(study)
    oreff = context.Output(reff_out)
    oieff = context.Output(ieff_out)
    for lane, rep, out in (("reff", ireff, oreff), ("ieff", iieff, oieff)):
        context.ExecWithContainer(
            image=env,
            cmd=f"""python {istud.container} run --lane {lane} \
                --report {rep.container} --nulls-dir {inul.container} \
                --faa {ifaa.container} --family-spec {ifam.container} \
                --out {out.container}""",
        )
    return ExecutionResult(
        manifest=[{reff_out: oreff.local, ieff_out: oieff.local}],
        success=oreff.local.exists() and oieff.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=2)),
)
