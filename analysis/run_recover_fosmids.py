#!/usr/bin/env python
"""Recover putative fosmid inserts from the pooled scadc reads -- WITHOUT reassembly.

Adaptive mid-start: the per-pool assemblies are pre-staged as sequences::assembly,
so the planner begins at clustering and never invokes an assembler.  DAG:

  recovery_experiment
    + sequences::assembly (per pool, spades_meta + megahit_default)
        -> cluster_contigs -> fosmids::clustered_contigs (+ cluster_membership)
    + sequences::clean_short_reads (one file per pool)
        + clustered_contigs -> pool_coverage -> fosmids::pool_coverage_profiles
  clustered_contigs + profiles -> chimera_split -> fosmids::split_contigs
  split_contigs     + profiles -> coverage_trim -> fosmids::putative_inserts (+ report)

A single fosmids::recovery_experiment node groups the whole run; every assembly and
pool-reads file is parented to it, so each transform (group_by=experiment) sees all
of them via InputGroup in one job.

Usage:
  python run_recover_fosmids.py               # generate + render + run (docker)
  python run_recover_fosmids.py --generate    # plan + render DAG only (no containers)

Run under the msm env (nextflow + docker + dot on PATH):
  /home/tony/lib/miniforge3/envs/msm/bin on PATH.
"""
import sys
import shutil
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

import pandas as pd

from metasmith.python_api import (
    Agent, ContainerRuntime, Source,
    DataInstanceLibrary, DataTypeLibrary, TransformInstanceLibrary,
    TargetBuilder, Resources, Size, Duration,
)

MLIB = Path(__file__).resolve().parent.parent
SCADC_DATA = Path("/home/tony/agentic_workspace/projects/scadc/recover-fosmids/.awm/data")
ASM_DIR = SCADC_DATA / "raw_runs" / "scadc_fabfos_assemblies"
KAT_DIR = SCADC_DATA / "raw_runs" / "kat_readmap"
POOL_LINEAGE = SCADC_DATA / "pool_lineage.csv"

# scadc used exactly these two assemblers for the 191-contig clustering
ASSEMBLERS = ["spades_meta", "megahit_default"]

STAGING = MLIB / ".awm" / "data" / "runs" / "recover_fosmids"
READS_DIR = STAGING / "pool_reads"

DOMAINS = [
    "assembly", "fosmids", "functionalAnnotation", "logistics",
    "metabolicModelling", "metagenomics", "pangenome",
    "responseSurface", "transcriptomics",
]

ASSEMBLER_TOKENS = ("megahit", "spades", "flye", "hifiasm", "miniasm", "assembl")


def materialize_pool_reads() -> dict[str, Path]:
    """Concatenate each pool's barcode reads (R1+R2) into one <pool>.fq.

    Pool 01 carries 3 barcodes; concatenating them sums their depth, matching
    scadc coverage_profile.build_profiles (which '+=' accumulates the barcodes).
    Reads are the immutable kat-unmapped fastqs; we only read them.
    """
    READS_DIR.mkdir(parents=True, exist_ok=True)
    lineage = pd.read_csv(POOL_LINEAGE, dtype=str)
    pool_fastqs: dict[str, Path] = {}
    for _, r in lineage.iterrows():
        pool = str(r["pool"]).strip()
        barcodes = [b.strip() for b in str(r["barcode"]).split(";") if b.strip()]
        out_fq = READS_DIR / f"{pool}.fq"
        pool_fastqs[pool] = out_fq
        if out_fq.exists() and out_fq.stat().st_size > 0:
            continue
        srcs = []
        for bc in barcodes:
            for mate in ("1", "2"):
                p = KAT_DIR / bc / f"{bc}_unmapped_{mate}.fq"
                if p.exists():
                    srcs.append(p)
                else:
                    print(f"  WARNING: missing reads {p}")
        if not srcs:
            print(f"  WARNING: pool {pool} has no reads; skipping")
            pool_fastqs.pop(pool, None)
            continue
        tmp = out_fq.with_suffix(".fq.partial")
        with open(tmp, "wb") as w:
            for p in srcs:
                with open(p, "rb") as rd:
                    shutil.copyfileobj(rd, w, length=16 * 1024 * 1024)
        tmp.rename(out_fq)
        print(f"  pool {pool}: {len(srcs)} files -> {out_fq.name}")
    return pool_fastqs


def build_inputs(materialize_reads: bool) -> DataInstanceLibrary:
    inputs_dir = STAGING / "inputs.xgdb"
    if inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    inputs = DataInstanceLibrary(inputs_dir)
    inputs.Purge()
    inputs.AddTypeLibrary(namespace="sequences",
                          lib=DataTypeLibrary.Load(MLIB / "data_types/sequences.yml"))
    inputs.AddTypeLibrary(namespace="fosmids",
                          lib=DataTypeLibrary.Load(MLIB / "data_types/fosmids.yml"))

    experiment = inputs.AddValue(
        "recovery_experiment.txt", "scadc_recover_fosmids",
        "fosmids::recovery_experiment",
    )

    n_asm = 0
    for asm_file in sorted(ASM_DIR.glob("*.fna")):
        stem = asm_file.stem
        if "__" not in stem:
            continue
        _pool, assembler = stem.split("__", 1)
        if assembler not in ASSEMBLERS:
            continue
        inputs.AddItem(asm_file, "sequences::assembly", parents={experiment})
        n_asm += 1
    assert n_asm > 0, f"no assemblies matched {ASSEMBLERS} under {ASM_DIR}"
    print(f"staged {n_asm} assemblies ({', '.join(ASSEMBLERS)})")

    if materialize_reads:
        pool_fastqs = materialize_pool_reads()
    else:
        # plan-only: reference intended paths without concatenating 3.9G
        lineage = pd.read_csv(POOL_LINEAGE, dtype=str)
        pool_fastqs = {str(r["pool"]).strip(): READS_DIR / f'{str(r["pool"]).strip()}.fq'
                       for _, r in lineage.iterrows()}

    for pool, fq in sorted(pool_fastqs.items()):
        inputs.AddItem(fq, "sequences::clean_short_reads", parents={experiment})
    print(f"staged {len(pool_fastqs)} pool read sets")

    inputs.Save()
    return inputs


def load_resources_and_transforms():
    resources = [DataInstanceLibrary.Load(MLIB / f"resources/{n}") for n in ["containers", "lib"]]
    transforms = [TransformInstanceLibrary.Load(MLIB / f"transforms/{d}") for d in DOMAINS]
    return resources, transforms


def assert_no_assembler(task):
    for step in task.plan.steps:
        name = Path(step.transform._path).stem.lower()
        if any(tok in name for tok in ("megahit", "spades", "flye", "hifiasm", "miniasm")):
            raise AssertionError(f"plan contains an assembler step [{name}] -- mid-start failed")


def main():
    generate_only = "--generate" in sys.argv or "--generate-only" in sys.argv

    print("=== Building input library ===")
    inputs = build_inputs(materialize_reads=not generate_only)

    print("\n=== Loading resources & transforms ===")
    resources, transforms = load_resources_and_transforms()

    print("\n=== Setting up local agent (docker) ===")
    home = Source.FromLocal(STAGING / "agent_home")
    agent = Agent(home=home, runtime=ContainerRuntime.DOCKER)

    print("\n=== Generating workflow ===")
    targets = TargetBuilder()
    targets.Add("fosmids::putative_inserts")
    targets.Add("fosmids::putative_insert_report")
    task = agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("fosmids::recovery_experiment")),
        resources=resources + [inputs],
        transforms=transforms,
        targets=targets,
    )
    if not task.ok:
        print(f"FAILED: {task}")
        sys.exit(1)
    print(f"Plan has {len(task.plan.steps)} steps")
    for step in task.plan.steps:
        name = Path(step.transform._path).stem
        prods = [i.dtype_name for g in step.produces for i in g]
        print(f"  Step {step.order}: {name} -> {prods}")
    assert_no_assembler(task)
    print("OK: no assembler step in plan (adaptive mid-start confirmed)")

    print("\n=== Rendering DAG ===")
    dag_dir = MLIB / "reports"
    dag_dir.mkdir(parents=True, exist_ok=True)
    try:
        task.plan.RenderDAG(dag_dir / "recover_fosmids_dag.svg")
        print(f"DAG -> {dag_dir / 'recover_fosmids_dag.svg'}")
    except Exception as e:
        print(f"DAG rendering failed (non-fatal): {e}")

    if generate_only:
        print("\n[--generate] plan + DAG only; skipping stage/run.")
        return

    print("\n=== Staging workflow ===")
    agent.Deploy()
    agent.StageWorkflow(task, on_exist="update", verify_external_paths=False)

    print("\n=== Running workflow (docker) ===")
    # The transforms declare generous memory (32 GB) for large datasets; the local
    # nextflow executor here caps at ~8 GB, so scale every step to fit this host.
    # The clustered set is modest (~hundreds of contigs) so 6 GB / 8 cpus is ample.
    agent.RunWorkflow(
        task,
        resource_overrides={
            "all": Resources(cpus=8, memory=Size.GB(6), duration=Duration(hours=6)),
        },
    )

    print("\n=== Waiting for completion ===")
    workspace = STAGING / "agent_home" / "runs" / task._key
    internals = workspace / "_metasmith"
    t0 = time.time()
    timeout = 6 * 3600
    while time.time() - t0 < timeout:
        if internals.exists():
            logs = sorted(p for p in internals.glob("logs.*") if "latest" not in p.name)
            if logs:
                main_log = logs[-1] / "main.log"
                if main_log.exists():
                    text = main_log.read_text(errors="ignore")
                    if "run completed at" in text:
                        break
                    if "ERROR" in text and "nextflow" in text.lower():
                        print("".join(text.splitlines(keepends=True)[-40:]))
                        sys.exit(1)
        time.sleep(10)
    else:
        print(f"TIMEOUT after {timeout}s")
        sys.exit(1)

    agent.CheckWorkflow(task)
    print(f"\nRun completed in {(time.time()-t0)/60:.1f}min")

    results_dir = workspace / "results"
    print(f"\n=== Results under {results_dir} ===")
    for bucket in ("fosmids-putative_inserts", "fosmids-putative_insert_report",
                   "fosmids-clustered_contigs", "fosmids-split_contigs"):
        d = results_dir / bucket
        if d.exists():
            for p in d.iterdir():
                print(f"  {bucket}: {p} ({p.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
