"""ECSPr significance: reff/ieff reports + staged null + fosmid proteins
-> reff_significance + ieff_significance.

Thin wrapper around `resources/lib/ecspr_significance.py run`. Scores each observed
per-fosmid delta against the staged null with the FULL-MIXTURE survival function,
size-matched by SF interpolation across the two flanking draw-size anchors. Runs both
lanes. The proteins FASTA supplies the per-fosmid ORF counts for size matching.

The scorer DERIVES its draw sizes from the staged `ecspr::frozen_null` directory and
holds no size constant of its own -- so the sizes are whatever the caller staged, and
the caller is responsible for staging a CURATED directory built from an explicit list.
Do not point this at a raw cache: caches accumulate retired sizes beside canonical
ones, and a glob would silently widen the basis.

Parity: `ecspr_significance.py selftest` reproduces the canonical scadc table on both
lanes to ~1e-16, including the split-contig ids that the retired `\\w`-based ORF regex
dropped silently. It is CPU-only and runs in seconds -- gate every commit with it.
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
reff     = model.AddRequirement(lib.GetType("ecspr::reff_axes_report"), parents={exp})
ieff     = model.AddRequirement(lib.GetType("ecspr::ieff_axes_report"), parents={exp})
proteins = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
nulls    = model.AddRequirement(lib.GetType("ecspr::frozen_null"))
env      = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
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
