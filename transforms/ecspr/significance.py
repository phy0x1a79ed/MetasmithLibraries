"""ECSPr significance: reff/ieff reports + reused frozen null + fosmid proteins
-> reff_significance + ieff_significance.

Thin wrapper around `resources/lib/ecspr_significance.py run` (selftest reproduces
scadc sig_emp_{reff,ieff}.tsv to 1e-16). Empirical-CCDF exceedance vs the REUSED
frozen metagenome null (never recomputed), BH-q within (element, null-style),
size-matched to nearest N in {21,34,56} by ORF count. Runs both lanes. The
proteins FASTA (fosmid ORFs) supplies the per-fosmid ORF counts for size-matching.
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
reff     = model.AddRequirement(lib.GetType("ecspr::reff_axes_report"), parents={exp})
ieff     = model.AddRequirement(lib.GetType("ecspr::ieff_axes_report"), parents={exp})
proteins = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
nulls    = model.AddRequirement(lib.GetType("ecspr::frozen_null"))
env      = model.AddRequirement(lib.GetType("envs::ecspr.condaenv"))
sig      = model.AddRequirement(lib.GetType("lib::ecspr_significance.py"))
reff_sig = model.AddProduct(lib.GetType("ecspr::reff_significance"))
ieff_sig = model.AddProduct(lib.GetType("ecspr::ieff_significance"))

def protocol(context: ExecutionContext):
    ireff = context.Input(reff)
    iieff = context.Input(ieff)
    ifaa  = context.Input(proteins)
    inul  = context.Input(nulls)
    isig  = context.Input(sig)
    oreff = context.Output(reff_sig)
    oieff = context.Output(ieff_sig)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {isig.container} run --lane reff \
            --report {ireff.container} --nulls-dir {inul.container} \
            --faa {ifaa.container} --out {oreff.container}""",
    )
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {isig.container} run --lane ieff \
            --report {iieff.container} --nulls-dir {inul.container} \
            --faa {ifaa.container} --out {oieff.container}""",
    )
    return ExecutionResult(
        manifest=[{reff_sig: oreff.local, ieff_sig: oieff.local}],
        success=oreff.local.exists() and oieff.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=2)),
)
