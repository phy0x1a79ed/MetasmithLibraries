#!/usr/bin/env python
"""FabFos analysis E2E -- fresh 4-lane annotation + ECSPr on the 132 reference inserts.

Mid-start driver: begins from the canonical 132 >=29kb reference inserts and runs the
whole analysis half through metasmith in **mamba runtime** (tools in conda envs, not
containers), **serially** (nextflow local executor, queueSize=1). DAG:

  recovery_experiment
    + fosmids::reference_inserts (the canonical 132 fasta)
        -> orfcall_inserts -> sequences::open_reading_frames (+ gff)
    open_reading_frames -> kofam_lane   -> kofam_hits   (GATED: KOfam HMM DB + fabfos-kofam)
                        -> dl_ec_lane   -> dlec_pred    (GATED: EZpred/ESM-C weights + fabfos-ml + GPU)
                        -> uniref_lane  -> uniref_hits  (GATED: UniRef50 .dmnd 24GB + fabfos-annot)
                        -> proteinbert_embed -> protein_embeddings (GATED: ProteinBERT + fabfos-ml)
    reac_prop  -> build_ec_bridge      -> ec_to_mnxr
    reac_xref + rhea2uniprot{,_trembl} -> build_uniprot_bridge -> uniprot_to_mnxr
    kofam_hits + dlec_pred + uniref_hits + {ko,ec,uniprot}_to_mnxr
        -> compile_evidence -> evidence_table
    evidence_table -> evidence_weights -> reaction_catalog / addition_weights
    evidence_weights + mnx_bipartite + biomass_axes -> base_graphs
    addition_weights + base_graphs + mnx_bipartite + biomass_axes -> solve -> reff/ieff reports
    reports + open_reading_frames + frozen_null -> significance -> reff/ieff_significance
    ... + addition_weights -> ablation -> ablation_importance

Reuse policy: ONLY the frozen metagenome null (ecspr::frozen_null) and the
reference tables (MetaNetX, Rhea, ko_to_mnxr bridge, mnx_bipartite, biomass_axes)
are reused. The fosmid annotations are produced FRESH by the gated annotator lanes
-- scadc's precomputed fosmid annotations are NOT reused for the deliverable.

Usage:
  python run_fabfos_e2e.py --generate   # plan + render the full DAG (no staging, no run)
  python run_fabfos_e2e.py              # stage + run serially in mamba mode (needs the GATES)

Run under an env with nextflow + dot on PATH (/home/tony/lib/miniforge3/envs/msm/bin),
and with the fabfos-* conda envs provisioned (mamba mode does NOT create them).
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

from metasmith.python_api import (
    Agent, Runtime, Source,
    DataInstanceLibrary, DataTypeLibrary, TransformInstanceLibrary,
    TargetBuilder, Resources, Size, Duration,
)

MLIB = Path(__file__).resolve().parent.parent
STAGING = MLIB / ".awm" / "data" / "runs" / "fabfos_e2e"

# ---- canonical fresh input ------------------------------------------------
CANON = Path("/home/tony/agentic_workspace/data/scadc/fabfos_2026")
REFERENCE_INSERTS = CANON / "putative_inserts_132_ge29kb.fna"

# ---- reused reference artifacts (staged inputs; never recomputed) ----------
D = Path("/home/tony/agentic_workspace/data/scadc")
MM = Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling"
          "/main/metabolic-modelling")
RN_CACHE = MM / "04_reaction_network" / "cache"
REFS = {
    # MNXR bridges / MetaNetX / Rhea reference tables
    "functional_annotation::ko_to_mnxr":      MM / "_reference_try1/betweenness/cache/ko_to_mnxr.tsv",
    "ecspr::metanetx_reac_prop":              D / "references/metanetx/reac_prop.tsv",
    "ecspr::metanetx_chem_prop":              D / "references/metanetx/chem_prop.tsv",
    "functional_annotation::metanetx_reac_xref": D / "references/metanetx/reac_xref.tsv",
    "functional_annotation::rhea2uniprot":       D / "references/rhea/rhea2uniprot.tsv",
    "functional_annotation::rhea2uniprot_trembl": D / "references/rhea/rhea2uniprot_trembl.tsv.gz",
    # ECSPr network reference artifacts (biomass axes = set2cat)
    "ecspr::biomass_axes":                    RN_CACHE / "biomass_dag_axes_set2cat.json",
}
# Directory-typed reused refs -- curated symlink dirs are built at stage time.
BIPARTITE_ELEMENTS = ["C", "N", "S", "P"]
FROZEN_NULL_FILES = [f"{lane}_null_canonical_N{n}.tsv"
                     for lane in ("reff", "ieff") for n in (21, 34, 56)]

# ---- GATED external annotator DBs / model weights -------------------------
# Provide these to run FRESH end-to-end. Intended default paths shown; override
# by editing here or symlinking into place. Planning (--generate) does not need
# them to exist; a real run does.
GATES = {
    "functional_annotation::kofam_profiles":   D / "references/kofam/profiles",
    "functional_annotation::kofam_ko_list":    D / "references/kofam/ko_list",
    "functional_annotation::uniref50_dmnd":    D / "references/uniref50/uniref50.dmnd",
    "functional_annotation::dlec_model":       D / "references/ezpred_esmc_600m",
    "functional_annotation::proteinbert_model": D / "references/proteinbert",
    "functional_annotation::reference_label_pool": D / "metabolic_modelling/embed_transfer_pool",
}

# Only the mamba-native fabfos analysis domains. Deliberately NOT loading
# metagenomics/assembly (their docker prodigal etc. would be picked structurally
# for open_reading_frames and break under Runtime.MAMBA, which reads a container
# requirement's file CONTENT as the conda env name).
DOMAINS = ["fosmids", "functionalAnnotation", "ecspr"]


def build_dir_ref(name: str, files: list[Path]) -> Path:
    """Curate a directory of symlinks for a DIRECTORY-typed reused ref."""
    d = STAGING / "refs" / name
    d.mkdir(parents=True, exist_ok=True)
    for f in files:
        link = d / f.name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(f)
    return d


def build_inputs(stage_files: bool) -> DataInstanceLibrary:
    import shutil
    inputs_dir = STAGING / "inputs.xgdb"
    if inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    inputs = DataInstanceLibrary(inputs_dir)
    inputs.Purge()
    for ns in ("sequences", "fosmids", "functional_annotation", "ecspr"):
        inputs.AddTypeLibrary(namespace=ns,
                              lib=DataTypeLibrary.Load(MLIB / f"data_types/{ns}.yml"))

    exp = inputs.AddValue("recovery_experiment.txt", "fabfos_2026_analysis",
                          "fosmids::recovery_experiment")

    # fresh input: the canonical 132 reference inserts
    inputs.AddItem(REFERENCE_INSERTS, "fosmids::reference_inserts", parents={exp})

    # reused single-file reference tables
    for type_name, path in REFS.items():
        inputs.AddItem(path, type_name)

    # reused directory refs (curated symlink dirs)
    bip_dir = build_dir_ref("mnx_bipartite",
                            [RN_CACHE / f"mnx_bipartite_{e}.pkl" for e in BIPARTITE_ELEMENTS])
    inputs.AddItem(bip_dir, "ecspr::mnx_bipartite")
    null_dir = build_dir_ref("frozen_null", [RN_CACHE / f for f in FROZEN_NULL_FILES])
    inputs.AddItem(null_dir, "ecspr::frozen_null")

    # gated annotator DBs / model weights
    for type_name, path in GATES.items():
        inputs.AddItem(path, type_name)

    inputs.Save()
    return inputs


def main():
    generate_only = "--generate" in sys.argv

    print("=== staging inputs ===")
    inputs = build_inputs(stage_files=not generate_only)

    print("=== loading resources & transforms ===")
    resources = [DataInstanceLibrary.Load(MLIB / f"resources/{n}")
                 for n in ("containers", "envs", "lib")]
    transforms = [TransformInstanceLibrary.Load(MLIB / f"transforms/{d}") for d in DOMAINS]

    print("=== agent (mamba runtime) ===")
    agent = Agent(home=Source.FromLocal(STAGING / "agent_home"), runtime=Runtime.MAMBA)

    print("=== generating workflow ===")
    targets = TargetBuilder()
    targets.Add("ecspr::reff_significance")
    targets.Add("ecspr::ieff_significance")
    targets.Add("ecspr::ablation_importance")
    targets.Add("ecspr::reaction_catalog")
    task = agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("fosmids::recovery_experiment")),
        resources=resources + [inputs],
        transforms=transforms,
        targets=targets,
    )
    if not task.ok:
        print(f"PLANNING FAILED: {task}")
        sys.exit(1)
    print(f"plan has {len(task.plan.steps)} steps")
    for step in task.plan.steps:
        name = Path(step.transform._path).stem
        prods = [i.dtype_name for g in step.produces for i in g]
        print(f"  step {step.order}: {name} -> {prods}")

    print("=== rendering DAG ===")
    dag = MLIB / "reports" / "fabfos_e2e_dag.svg"
    dag.parent.mkdir(parents=True, exist_ok=True)
    try:
        task.plan.RenderDAG(dag)
        print(f"DAG -> {dag}")
    except Exception as e:
        print(f"DAG render failed (non-fatal): {e}")

    if generate_only:
        print("[--generate] plan + DAG only; skipping stage/run.")
        return

    print("=== staging workflow ===")
    agent.Deploy()
    agent.StageWorkflow(task, on_exist="update", verify_external_paths=False)

    print("=== running (mamba, serial queueSize=1) ===")
    agent.RunWorkflow(
        task,
        params=dict(executor=dict(queueSize=1), process=dict(maxForks=1)),
        resource_overrides={
            "all": Resources(cpus=8, memory=Size.GB(6), duration=Duration(hours=12)),
        },
    )
    print("submitted; monitor at",
          STAGING / "agent_home" / "runs" / task._key / "_metasmith" / "logs.latest")


if __name__ == "__main__":
    main()
