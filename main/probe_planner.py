#!/usr/bin/env python3
"""Plan-only probe: run GenerateWorkflow with various target sets and report
which transforms got selected. No staging, no remote work.

Designed to isolate why `--only all` adds `downloadESMFoldWeights` but
`--only esmc` doesn't, despite both having pre-staged weight tarballs.
"""
import sys
from pathlib import Path

MSM_SRC = Path("/home/tony/agentic_workspace/projects/metasmith/dev/src")
sys.path.insert(0, str(MSM_SRC))

from metasmith.python_api import (
    Agent, SshSource, DataInstanceLibrary, TransformInstanceLibrary,
    TargetBuilder, ContainerRuntime,
)

ROOT = Path(__file__).resolve().parent
DL_LIB = ROOT.parent
CACHE_DIR = ROOT / "cache"

HPC_HOST = "fir"
HPC_MSM_HOME = Path("/scratch/phyberos/metasmith")

# Targets identical to launch_dl_embeddings.py TARGETS table.
TARGETS = {
    "esmc":         "annotation::esm_c_embeddings",
    "ankh":         "annotation::ankh_embeddings",
    "prott5":       "annotation::prott5_embeddings",
    "esmfold":      "sequences::predicted_structures",
    "foldseek_3di": "sequences::structure_3di_tokens",
    "saprot":       "annotation::saprot_embeddings",
}

WEIGHTS = {
    "esmc":   ("/scratch/phyberos/dl_testing_claude/weights/esmc_300m.tgz",  "ref::esm_c_300m_weights"),
    "ankh":   ("/scratch/phyberos/dl_testing_claude/weights/ankh_base.tgz",  "ref::ankh_base_weights"),
    "prott5": ("/scratch/phyberos/dl_testing_claude/weights/prott5_xl.tgz",  "ref::prott5_xl_uniref50_weights"),
    "saprot": ("/scratch/phyberos/dl_testing_claude/weights/saprot_650m.tgz","ref::saprot_650m_weights"),
    "esmfold":("/scratch/phyberos/dl_testing_claude/weights/esmfold_v1.tgz", "ref::esmfold_weights"),
}

ORFS_PATH = Path("/scratch/phyberos/dl_testing_claude/inputs/scadc_metag.faa")


def get_agent():
    return Agent(
        home=SshSource(host=HPC_HOST, path=HPC_MSM_HOME).AsSource(),
        runtime=ContainerRuntime.APPTAINER,
        setup_commands=["module load apptainer"],
    )


def build_library(weights_keys: list[str], lib_name: str) -> DataInstanceLibrary:
    """Build a fresh local .xgdb with ORFs + listed weights."""
    in_dir = CACHE_DIR / f"probe_{lib_name}.xgdb"
    if in_dir.exists():
        import shutil
        shutil.rmtree(in_dir)
    lib = DataInstanceLibrary(in_dir)
    lib.Purge()
    lib.AddTypeLibrary(DL_LIB / "data_types" / "sequences.yml")
    lib.AddTypeLibrary(DL_LIB / "data_types" / "annotation.yml")
    lib.AddTypeLibrary(DL_LIB / "data_types" / "ref.yml")
    lib.AddItem(ORFS_PATH, "sequences::orfs")
    for w in weights_keys:
        path, type_str = WEIGHTS[w]
        lib.AddItem(Path(path), type_str)
    lib.Save()
    return lib


def probe(label: str, target_keys: list[str], weight_keys: list[str],
          max_iter: int = 1024, max_refine: int = 256, seed: int = 42):
    """Plan one configuration; report which transform names appear in the DAG."""
    print(f"\n{'='*70}\n[{label}] targets={target_keys}  weights={weight_keys}  "
          f"iter={max_iter} refine={max_refine} seed={seed}\n{'='*70}")

    inputs = build_library(weight_keys, label)
    containers = DataInstanceLibrary.Load(DL_LIB / "resources" / "containers")
    transforms = [
        TransformInstanceLibrary.Load(DL_LIB / "transforms" / "functionalAnnotation"),
        TransformInstanceLibrary.Load(DL_LIB / "transforms" / "logistics"),
    ]

    tb = TargetBuilder()
    for t in target_keys:
        tb.Add(TARGETS[t])

    smith = get_agent()
    task = smith.GenerateWorkflow(
        samples=list(inputs.AsSamples("sequences::orfs")),
        resources=[containers, inputs],
        transforms=transforms,
        targets=tb,
        max_iter=max_iter, max_refine=max_refine, seed=seed,
    )
    if not task.ok or len(task.plan.steps) == 0:
        print(f"  PLAN FAILED ok={task.ok} steps={len(task.plan.steps)}")
        return

    # Build Transform.key → TransformInstance.name map
    tr_name_by_key: dict[str, str] = {}
    for trlib in transforms:
        for path, ti in trlib.IterateTransforms():
            tr_name_by_key[ti.model.key] = ti.name or str(path)

    print(f"  OK steps={len(task.plan.steps)} key={task.GetKey()}")
    for i, step in enumerate(task.plan.steps, 1):
        ti = step.transform   # TransformInstance, not Transform
        key = ti.model.key if hasattr(ti, "model") else getattr(ti, "_key", "?")
        name = tr_name_by_key.get(key, getattr(ti, "name", None) or f"<unmapped key={key}>")
        print(f"    {i:2d}. {name}")


if __name__ == "__main__":
    case = sys.argv[1] if len(sys.argv) > 1 else "esmfold_only"

    if case == "esmfold_only":
        probe("esmfold_only", ["esmfold"], ["esmfold"])
    elif case == "saprot_only":
        # saprot transitively needs esmfold→foldseek_3di chain + its own weights + esmfold weights
        probe("saprot_only", ["saprot"], ["saprot", "esmfold"])
    elif case == "esmc_baseline":
        probe("esmc_baseline", ["esmc"], ["esmc"])
    elif case == "all":
        probe("all", list(TARGETS.keys()), ["esmc", "ankh", "prott5", "saprot", "esmfold"])
    elif case == "esmc_plus_esmfold":
        # does adding esmfold to esmc cause download to appear for esmfold but not esmc?
        probe("esmc_plus_esmfold", ["esmc", "esmfold"], ["esmc", "esmfold"])
    elif case == "all_minus_saprot":
        # everything except saprot/foldseek_3di — strips the structural chain
        probe("all_minus_saprot", ["esmc", "ankh", "prott5", "esmfold"],
              ["esmc", "ankh", "prott5", "esmfold"])
    elif case == "all_no_foldseek_target":
        # saprot transitively pulls foldseek_3di; drop it as explicit target
        probe("all_no_foldseek_target",
              ["esmc", "ankh", "prott5", "esmfold", "saprot"],
              ["esmc", "ankh", "prott5", "saprot", "esmfold"])
    elif case == "esmfold_plus_foldseek":
        # minimal: do explicit foldseek_3di target trigger download?
        probe("esmfold_plus_foldseek", ["esmfold", "foldseek_3di"], ["esmfold"])
    elif case == "all_seed_variants":
        # same as 'all' but vary nothing — sanity that seed=42 is stable
        probe("all_seed_stable",
              ["esmc", "ankh", "prott5", "esmfold", "foldseek_3di", "saprot"],
              ["esmc", "ankh", "prott5", "saprot", "esmfold"])
    elif case == "all_no_esmfold_target":
        # saprot transitively needs esmfold; drop esmfold as explicit target
        probe("all_no_esmfold_target",
              ["esmc", "ankh", "prott5", "foldseek_3di", "saprot"],
              ["esmc", "ankh", "prott5", "saprot", "esmfold"])
    elif case == "esmfold_plus_saprot":
        # minimal hypothesized trigger: just esmfold + saprot
        probe("esmfold_plus_saprot", ["esmfold", "saprot"], ["esmfold", "saprot"])
    elif case == "esmfold_plus_saprot_noweights":
        probe("esmfold_plus_saprot_noweights", ["esmfold", "saprot"], ["esmfold"])
    elif case == "1enc_struct":
        probe("1enc_struct", ["esmc", "esmfold", "saprot"], ["esmc", "esmfold", "saprot"])
    elif case == "2enc_struct":
        probe("2enc_struct", ["esmc", "ankh", "esmfold", "saprot"],
              ["esmc", "ankh", "esmfold", "saprot"])
    elif case == "3enc_struct":
        probe("3enc_struct", ["esmc", "ankh", "prott5", "esmfold", "saprot"],
              ["esmc", "ankh", "prott5", "esmfold", "saprot"])
    elif case == "3enc_no_saprot":
        probe("3enc_no_saprot", ["esmc", "ankh", "prott5", "esmfold"],
              ["esmc", "ankh", "prott5", "esmfold"])
    elif case == "all_big_iter":
        # same as 'all' but 8x MCTS budget
        probe("all_big_iter",
              ["esmc", "ankh", "prott5", "esmfold", "foldseek_3di", "saprot"],
              ["esmc", "ankh", "prott5", "saprot", "esmfold"],
              max_iter=8192, max_refine=2048)
    elif case == "all_seeds":
        # vary seed to see if different rolls find the no-download plan
        for s in [1, 7, 13, 42, 99, 2024]:
            probe(f"all_seed{s}",
                  ["esmc", "ankh", "prott5", "esmfold", "foldseek_3di", "saprot"],
                  ["esmc", "ankh", "prott5", "saprot", "esmfold"],
                  seed=s)
    else:
        print(f"unknown case: {case}")
        sys.exit(1)
