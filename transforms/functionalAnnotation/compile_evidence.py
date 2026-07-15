"""Evidence compiler -- fold the annotator lanes + MNXR bridges -> evidence_table.

Thin wrapper around `resources/lib/fabfos_evidence.py compile`. Runs the canonical
lanes (kofam KO->MNXR, CLEAN EC->MNXR, uniref50 UniProt->MNXR) through their readers
and concatenates into the unified 8-col evidence table ECSPr consumes. Ported from
scadc 11_build_evidence_table_dlec.py, with CLEAN as the EC channel (replacing
EZpred/dl_ec; dl_ec is retained in the lib as an optional comparison lane via
`--dl-ec`). The lane-4 embed-transfer candidates are produced separately; fold them
in via the CLI's --embed once staged.

ko_to_mnxr is a reused reference table (staged input); ec_to_mnxr and uniprot_to_mnxr
are built fresh by the bridge transforms.

WHY `source` IS A STAGED INPUT AND NOT A PARAMETER
--------------------------------------------------
It used to be `context.params.get("source", "fosmid")`, which is two silent bugs
wearing one coat. params are populated ONLY from the runtime resources line
(cpus/memory/attempt), so the value was ALWAYS "fosmid" no matter what a caller
passed -- and params never enter the task hash, so a second compile of a DIFFERENT
organism would have collided with the fosmid compile on one cache entry and been
served the fosmid table. `source` is the column every downstream consumer partitions
on, so those two failures compose into "this organism's reactions are now that
organism's reactions", with nothing raising. As a staged, content-hashed input the
tag enters the cache key: two sources cannot collide, and a missing tag is an error
rather than a default. This transform is single-source per invocation BY DESIGN --
the host and the inserts are separate compiles, joined downstream.
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
env      = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
src      = model.AddRequirement(lib.GetType("functional_annotation::evidence_source"),
                                parents={exp})
ev       = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
table    = model.AddProduct(lib.GetType("functional_annotation::evidence_table"))


def read_source(path):
    """`source: <tag>` -> the tag that labels every row of this compile.

    Fails loudly on an empty or malformed file. The value is the `source` column of
    the evidence table, and every consumer partitions on it -- so a wrong tag is not
    a cosmetic mislabel, it silently reassigns one organism's reactions to another.
    """
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        k, _, v = line.partition(":")
        if k.strip() == "source" and v.strip():
            return v.strip()
    raise AssertionError(f"evidence_source: no `source:` line in {path}")

def protocol(context: ExecutionContext):
    ikof = context.Input(kofam)
    icln = context.Input(clean)
    iuni = context.Input(uniref)
    ikb  = context.Input(ko_br)
    ieb  = context.Input(ec_br)
    iub  = context.Input(up_br)
    iev  = context.Input(ev)
    isrc = context.Input(src)
    otab = context.Output(table)
    source = read_source(isrc.local)
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
