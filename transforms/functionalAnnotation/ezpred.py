"""EZpred (DL-only fork) EC-number prediction transform -> ezpred_predictions.

Runs the EZpred deep-learning EC heads (kad-ecoli/EZpred, MIT) on ESM-C 600M
embeddings. EZpred's own path re-embeds each FASTA with its internal fasta2plm.py;
here we feed it the embeddings the `esm_c` transform already produced, so the
dataflow is esm_c -> ezpred and the expensive 600M pass runs once for both the
ESM-C GPR lane and EZpred (the original stub's TODO #1). Retyped/rewired from the
deep-learning functionalAnnotation/ezpred.py.

Inputs:
  - annotation::esm_c_embeddings : ESM-C 600M embedding vectors (from esm_c)
  - annotation::esm_c_index      : sequence_id -> row index for the embeddings
  - ref::ezpred_model            : the EZpred model bundle -- vendored EZpred tree
                                   (predict.py + the enzyme / non-enzyme MLP
                                   ensembles). It no longer needs the ESM-C
                                   weights, since embeddings are supplied.
  - env::ezpred.env              : the ESM-C image (esm SDK + torch + numpy/pandas
                                   /sklearn); EZpred has no image of its own. The
                                   EZpred-only dep (iterative-stratification) is
                                   pip-installed at start.

Output:
  - annotation::ezpred_predictions : CSV with columns
        sequence_id, ec_number, score, head_kind
    head_kind in {enzyme, nonenzyme} -- concat of DL1.tsv + DL2.tsv from
    predict.py --dl-only, melted to long form. The GPR mapper keeps the enzyme
    head, level-4 ECs above the DL_EC_SCORE_FLOOR.

RUNTIME CAVEAT (deferred): EZpred's heads were trained on the last-3-layer ESM-C
means, whereas esm_c emits a single mean-pooled vector. Feeding EZpred these
embeddings therefore needs either esm_c to also emit the per-layer means or a
predict.py patch to accept the single-vector feature -- a make-it-run task, out of
scope for the planning DAG. GATED: needs the ESM-C image + the ezpred_model bundle
staged; not run here.
"""
import os
from metasmith.python_api import *
from pathlib import Path

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image    = model.AddRequirement(lib.GetType("env::ezpred.env"))
ezmodel  = model.AddRequirement(lib.GetType("ref::ezpred_model"))
emb      = model.AddRequirement(lib.GetType("annotation::esm_c_embeddings"))
idx      = model.AddRequirement(lib.GetType("annotation::esm_c_index"))
out_pred = model.AddProduct(lib.GetType("annotation::ezpred_predictions"))


# Module-level constants -- same convention as esm_c.py.
TOP_TERMS = 500   # passed to predict.py --top_terms
# EZpred-only Python dep not in the ESM-C image (torch/numpy/pandas/sklearn/esm
# already live there); installed at protocol start.
EZPRED_PIP = "iterative-stratification==0.1.7"


WRAPPER = r'''
import os, sys, glob, subprocess, time, pickle
import numpy as np
import pandas as pd

# /work/EZpred is the bind-mounted model bundle (vendored EZpred tree + MLP heads).
EZPRED_DIR = "/work/EZpred"
EMB_PARQUET = os.environ["EZPRED_EMB"]
EMB_INDEX   = os.environ["EZPRED_INDEX"]
OUT_CSV     = os.environ["EZPRED_OUT_CSV"]
TOP_TERMS   = int(os.environ.get("EZPRED_TOP_TERMS", "500"))

# Build EZpred's precomputed-feature cache from the esm_c embeddings, so predict.py
# --dl-only skips its internal fasta2plm re-embed. EZpred consumes features keyed
# by sequence_id in the {name, mean: {layer: vec}, homo: {layer: vec}} shape.
idx = pd.read_csv(EMB_INDEX)
vecs = pd.read_parquet(EMB_PARQUET).to_numpy(dtype=np.float32)
feats = {}
for _, r in idx.iterrows():
    sid = str(r["sequence_id"])
    v = vecs[int(r["index"])]
    feats[sid] = {"name": sid, "mean": {0: v}, "homo": {0: v}}

work = os.path.join(EZPRED_DIR, "Data", "embeddings", "current_run")
os.makedirs(work, exist_ok=True)
feat_pkl = os.path.join(work, "precomputed_features.pkl")
with open(feat_pkl, "wb") as fh:
    pickle.dump(feats, fh)
print(f"staged {len(feats)} precomputed ESM-C features -> {feat_pkl}", flush=True)

t0 = time.time()
cmd = [
    "python", os.path.join(EZPRED_DIR, "predict.py"),
    "-w", work,
    "--features", feat_pkl,          # consume precomputed embeddings; bypass fasta2plm
    "--top_terms", str(TOP_TERMS),
    "--dl-only",
]
print("running:", " ".join(cmd), flush=True)
r = subprocess.run(cmd, cwd=EZPRED_DIR)
if r.returncode != 0:
    sys.exit(f"EZpred predict.py exited {r.returncode}")
print(f"predict.py finished in {time.time()-t0:.1f}s", flush=True)

# Concatenate DL1.tsv (enzyme) + DL2.tsv (nonenzyme) into a long CSV.
# predict.py writes rows = sequences, columns = EC terms; melt to long form.
out_rows = []
for kind, tsv_name in [("enzyme", "DL1.tsv"), ("nonenzyme", "DL2.tsv")]:
    path = os.path.join(work, tsv_name)
    if not os.path.exists(path):
        sys.exit(f"missing expected output: {path}")
    df = pd.read_csv(path, sep="\t")
    id_col = df.columns[0]
    long = df.melt(id_vars=[id_col], var_name="ec_number", value_name="score")
    long = long.rename(columns={id_col: "sequence_id"})
    long["head_kind"] = kind
    out_rows.append(long[["sequence_id", "ec_number", "score", "head_kind"]])
out = pd.concat(out_rows, ignore_index=True)
out.to_csv(OUT_CSV, index=False)
print(f"wrote {len(out)} rows to {OUT_CSV}", flush=True)
'''


def protocol(context: ExecutionContext):
    iemb  = context.Input(emb)
    iidx  = context.Input(idx)
    iez   = context.Input(ezmodel)
    ipred = context.Output(out_pred)

    wrapper = Path("ezpred_wrapper.py")
    with open(wrapper, "w") as f:
        f.write(WRAPPER)

    # iez.local is the model-bundle directory. Bind-mount it at /work/EZpred so
    # settings.py's root_dir resolves correctly.
    context.ExecWithEnv().ifContainerDo(
        env=image,
        binds=[
            (context.external_cwd/wrapper.name, f"/work/{wrapper.name}"),
            (iez.local, "/work/EZpred"),
        ],
        args=["--nv", "--env", f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES','')}"],
        cmd=f"""
            pip install --quiet --no-cache-dir {EZPRED_PIP}
            EZPRED_EMB={iemb.container} \
            EZPRED_INDEX={iidx.container} \
            EZPRED_OUT_CSV={ipred.container} \
            EZPRED_TOP_TERMS={TOP_TERMS} \
            python /work/{wrapper.name}
        """,
    )

    return ExecutionResult(
        manifest=[{out_pred: ipred.local}],
        success=ipred.local.exists(),
    )


# Just the EC-head MLP inference now (5 ensembles x 2 heads = 10 MLPs) -- the
# heavy ESM-C 600M embed pass is done upstream in esm_c and reused here.
TransformInstance(
    protocol=protocol,
    model=model,
    group_by=emb,
    resources=Resources(
        cpus=4,
        memory=Size.GB(16),
        duration=Duration(hours=2),
    ),
)
