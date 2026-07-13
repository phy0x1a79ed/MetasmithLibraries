"""Lane 2 (canonical) -- CLEAN contrastive EC prediction for the fosmid ORFs -> clean_pred.

CLEAN (Yu et al., "Enzyme function prediction using contrastive learning",
Science 2023) places a query enzyme by *function*: its supervised-contrastive
objective pulls convergent-function / divergent-sequence enzymes together, so it
reaches enzymes the homology channels (kofam KO, UniRef DIAMOND) miss. It is
sequence-input but homology-defeating by construction -- the reason it replaces
EZpred/`dl_ec` as the canonical EC lane for ECSPr.

Runs the baked `external_clean` image (CUDA torch + ESM-1b weights + CLEAN
maxsep assets under /app). ESM-1b (650M) mean-embeds the ORFs on the GPU, then
CLEAN max-separation inference (`CLEAN_infer_fasta.py`) compares each query to the
EC-cluster-center embeddings and emits, per ORF, the selected EC set with a
GMM-calibrated confidence per call. Output is standardized to the 3-col TSV the
evidence compiler's `read_clean` consumes:

    Query ID <TAB> Predicted EC number <TAB> clean_score

`clean_score` is CLEAN's raw (GMM-calibrated) maxsep value. CLEAN never abstains,
so `read_clean` gates the lane at the F1-optimal clean_score >= 0.01 floor
downstream -- do not threshold here. CLEAN writes relative to CWD and /app is
read-only under apptainer, so we run from a writable, bind-mounted /clean_ws that
symlinks the baked read-only assets.

GATED: needs the `clean.oci` image (~ships ESM-1b weights + pretrained bundle)
and a GPU to be practical.
"""
from pathlib import Path
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
image = model.AddRequirement(lib.GetType("containers::clean.oci"))
exp   = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
orfs  = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
pred  = model.AddProduct(lib.GetType("functional_annotation::clean_pred"))

# In-container wrapper: builds a writable workspace that symlinks the baked CLEAN
# assets, runs maxsep inference over the ORFs, and reshapes CLEAN's ragged maxsep
# CSV into the standardized 3-col TSV. Runs entirely offline (ESM-1b weights +
# pretrained bundle baked into the image).
WRAPPER = r'''
import argparse, os, subprocess, sys, shutil

p = argparse.ArgumentParser()
p.add_argument("--fasta", required=True, help="input ORF FASTA (container path)")
p.add_argument("--out", required=True, help="output standardized TSV (container path)")
p.add_argument("--app", default="/app", help="baked CLEAN install root")
p.add_argument("--workdir", default="/clean_ws", help="writable working dir (bind-mounted)")
a = p.parse_args()

APP, WORK, NAME = a.app, a.workdir, "orfs"  # CLEAN keys all relative paths off this basename

os.makedirs(os.path.join(WORK, "data", "inputs"), exist_ok=True)
os.makedirs(os.path.join(WORK, "data", "esm_data"), exist_ok=True)
os.makedirs(os.path.join(WORK, "results", "inputs"), exist_ok=True)

def link(src, dst):
    if not os.path.lexists(dst):
        os.symlink(src, dst)

# CLEAN reads ./data/pretrained/{split100.pth,100.pt,gmm_ensumble.pkl},
# ./data/split100.csv (EC label map), and runs ./esm/scripts/extract.py -- all
# relative to CWD.
link(os.path.join(APP, "data", "pretrained"), os.path.join(WORK, "data", "pretrained"))
link(os.path.join(APP, "data", "split100.csv"), os.path.join(WORK, "data", "split100.csv"))
link(os.path.join(APP, "esm"), os.path.join(WORK, "esm"))
shutil.copy(a.fasta, os.path.join(WORK, "data", "inputs", f"{NAME}.fasta"))

env = dict(os.environ)
env.setdefault("TORCH_HOME", "/opt/torch_cache")
print(f"[clean] running maxsep on {a.fasta} (cwd={WORK})", flush=True)
r = subprocess.run(
    [sys.executable, os.path.join(APP, "CLEAN_infer_fasta.py"), "--fasta_data", NAME],
    cwd=WORK, env=env,
)
if r.returncode != 0:
    sys.exit(f"[clean] CLEAN_infer_fasta.py failed (rc={r.returncode})")

# Parse CLEAN's ragged maxsep CSV -> standardized 3-col TSV.
# Format (no header): <seq_id>,EC:<ec>/<score>,EC:<ec>/<score>,...
res_csv = os.path.join(WORK, "results", "inputs", f"{NAME}_maxsep.csv")
if not os.path.exists(res_csv):
    sys.exit(f"[clean] expected maxsep output missing: {res_csv}")

n_rows = n_orfs = 0
with open(res_csv) as fin, open(a.out, "w") as fout:
    fout.write("Query ID\tPredicted EC number\tclean_score\n")
    for line in fin:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split(",")
        # take the first whitespace token of the FASTA header so the ORF id
        # matches the kofam/uniref lanes (which key on `>`-header field 0).
        sid = parts[0].strip().split()[0] if parts[0].strip() else parts[0].strip()
        n_orfs += 1
        for tok in parts[1:]:
            tok = tok.strip()
            if not tok:
                continue
            body = tok[3:] if tok.startswith("EC:") else tok   # "EC:3.6.1.43/8.06" -> "3.6.1.43/8.06"
            ec, score = (body.rsplit("/", 1) if "/" in body else (body, "NA"))
            fout.write(f"{sid}\t{ec.strip()}\t{score.strip()}\n")
            n_rows += 1

print(f"[clean] wrote {n_rows} EC calls across {n_orfs} ORFs -> {a.out}", flush=True)
'''


def protocol(context: ExecutionContext):
    import os
    iorfs = context.Input(orfs)
    opred = context.Output(pred)

    # writable workspace bound into the container (CLEAN writes relative to CWD;
    # /app is read-only under apptainer).
    context.LocalShell("mkdir -p clean_ws")
    script = Path("run_clean.py")
    with open(script, "w") as f:
        f.write(WRAPPER)

    context.ExecWithContainer(
        image=image,
        binds=[
            (context.external_cwd / "clean_ws", "/clean_ws"),
            (context.external_cwd / script.name, f"/work/{script.name}"),
        ],
        args=[
            "--nv",
            "--env", f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES','')}",
            "--env", "TORCH_HOME=/opt/torch_cache",
        ],
        cmd=f"""
            python /work/{script.name} \
                --fasta {iorfs.container} \
                --out {opred.container} \
                --app /app \
                --workdir /clean_ws
        """,
    )
    return ExecutionResult(
        manifest=[{pred: opred.local}],
        success=opred.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    # ESM-1b 650M dominates; embeds a fosmid-scale ORF set in minutes on a V100.
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=4)),
)
