"""Lane 2 (canonical) -- deep-learning EC prediction for the fosmid ORFs -> dlec_pred.

Fresh EZpred / ESM-C 600M EC-head inference (mamba, fabfos-ml env, GPU) over the
reference-insert ORFs. Emits a parquet with columns (sequence_id, ec_number, score,
head_kind) that `resources/lib/fabfos_evidence.py::read_dl_ec` filters to enzyme-head
level-4 ECs clearing DL_EC_SCORE_FLOOR (0.3); the compiler projects EC -> MNXR. This
is the canonical EC channel (deepEC is an optional comparison lane, not authored here).

GATED: needs the EZpred / ESM-C 600M EC-head weights (dlec_model dir) + the fabfos-ml
env (torch + CUDA) and a GPU to be practical. The `ezpred` CLI must be on PATH in
that env (packaged research tool -- provision alongside the weights).
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
exp   = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
orfs  = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
mdl   = model.AddRequirement(lib.GetType("functional_annotation::dlec_model"))
env   = model.AddRequirement(lib.GetType("envs::ml.condaenv"))
pred  = model.AddProduct(lib.GetType("functional_annotation::dlec_pred"))

def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    imdl  = context.Input(mdl)
    opred = context.Output(pred)
    device = context.params.get("device", "cuda")
    context.ExecWithContainer(
        image=env,
        cmd=f"""ezpred predict \
            --model {imdl.container} \
            --fasta {iorfs.container} \
            --device {device} \
            --out {opred.container}""",
    )
    return ExecutionResult(
        manifest=[{pred: opred.local}],
        success=opred.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=6)),
)
