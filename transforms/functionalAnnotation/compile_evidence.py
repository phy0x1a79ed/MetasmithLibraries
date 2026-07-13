"""Evidence compiler -- fold the annotator lanes + MNXR bridges -> evidence_table.

Thin wrapper around `resources/lib/fabfos_evidence.py compile`. Runs the canonical
lanes (kofam KO->MNXR, CLEAN EC->MNXR, uniref50 UniProt->MNXR) through their readers
and concatenates into the unified 8-col evidence table ECSPr consumes. Ported from
scadc 11_build_evidence_table_dlec.py, single-source (fosmid), with CLEAN as the EC
channel (replacing EZpred/dl_ec; dl_ec is retained in the lib as an optional
comparison lane via `--dl-ec`). The lane-4 embed-transfer candidates are produced
separately; fold them in via the CLI's --embed once staged.

ko_to_mnxr is a reused reference table (staged input); ec_to_mnxr and uniprot_to_mnxr
are built fresh by the bridge transforms.
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
kofam    = model.AddRequirement(lib.GetType("functional_annotation::kofam_hits"), parents={exp})
clean    = model.AddRequirement(lib.GetType("functional_annotation::clean_pred"), parents={exp})
uniref   = model.AddRequirement(lib.GetType("functional_annotation::uniref_hits"), parents={exp})
ko_br    = model.AddRequirement(lib.GetType("functional_annotation::ko_to_mnxr"))
ec_br    = model.AddRequirement(lib.GetType("functional_annotation::ec_to_mnxr"))
up_br    = model.AddRequirement(lib.GetType("functional_annotation::uniprot_to_mnxr"))
env      = model.AddRequirement(lib.GetType("envs::ecspr.condaenv"))
ev       = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
table    = model.AddProduct(lib.GetType("functional_annotation::evidence_table"))

def protocol(context: ExecutionContext):
    ikof = context.Input(kofam)
    icln = context.Input(clean)
    iuni = context.Input(uniref)
    ikb  = context.Input(ko_br)
    ieb  = context.Input(ec_br)
    iub  = context.Input(up_br)
    iev  = context.Input(ev)
    otab = context.Output(table)
    source = context.params.get("source", "fosmid")
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {iev.container} compile \
            --source {source} \
            --kofam {ikof.container} \
            --clean {icln.container} \
            --uniref50 {iuni.container} \
            --ko-to-mnxr {ikb.container} \
            --ec-to-mnxr {ieb.container} \
            --uniprot-to-mnxr {iub.container} \
            --out {otab.container}""",
    )
    return ExecutionResult(
        manifest=[{table: otab.local}],
        success=otab.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(minutes=30)),
)
