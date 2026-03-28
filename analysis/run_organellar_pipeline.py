#!/usr/bin/env python
"""Run the organellar transcriptomics pipeline on Sockeye via SLURM.

Pipeline:
  GenBank → genbank_to_reference → organellar_reference + organellar_gff
  organellar_reference + read_pairs → minimap2_rnaseq_align → organellar_bam
  organellar_bam + organellar_gff + experiment → organellar_count_matrix → organellar_gene_count_table
  organellar_gene_count_table + experiment → deseq2 → deseq2_results
  organellar_reference + assembly → blast_contamination_check → contamination_report

Targets: organellar_gene_count_table, deseq2_results, contamination_report
"""
import sys
import time
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

from pathlib import Path
from metasmith.python_api import (
    Agent, ContainerRuntime, Source, SshSource,
    DataInstanceLibrary, TransformInstanceLibrary,
    TargetBuilder, Resources, Size,
)

MLIB = Path(__file__).resolve().parent.parent
SOCKEYE_RAW = Path("/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/porphyridium/Raw-transcriptome-data")
SOCKEYE_ASSEMBLY = Path("/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/eguEpdhP-intermediates/assembly/1-1-1.f1CMorcneUoLGMna-O4PhHAkd.fna")
SOCKEYE_BAM_CACHE = Path("/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/organellar_bams")

# GenBank files for organellar genomes (on Sockeye arc storage)
SOCKEYE_GENBANK = Path("/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/organellar_genbank")

SAMPLES = {
    "POR-0-1": ("POR-0-1-090325_S59_L001_R1_001.fastq.gz", "POR-0-1-090325_S59_L001_R2_001.fastq.gz"),
    "POR-0-2": ("POR-0-2-090325_S60_L001_R1_001.fastq.gz", "POR-0-2-090325_S60_L001_R2_001.fastq.gz"),
    "POR-0-3": ("POR-0-3-090325_S61_L001_R1_001.fastq.gz", "POR-0-3-090325_S61_L001_R2_001.fastq.gz"),
    "POR-S-1": ("POR-S-1-090325_S62_L001_R1_001.fastq.gz", "POR-S-1-090325_S62_L001_R2_001.fastq.gz"),
    "POR-S-2": ("POR-S-2-090325_S67_L001_R1_001.fastq.gz", "POR-S-2-090325_S67_L001_R2_001.fastq.gz"),
    "POR-S-3": ("POR-S-3-090325_S68_L001_R1_001.fastq.gz", "POR-S-3-090325_S68_L001_R2_001.fastq.gz"),
    "POR-8-1": ("POR-8-1-090325_S66_L001_R1_001.fastq.gz", "POR-8-1-090325_S66_L001_R2_001.fastq.gz"),
    "POR-8-2": ("POR-8-2-090325_S69_L001_R1_001.fastq.gz", "POR-8-2-090325_S69_L001_R2_001.fastq.gz"),
    "POR-8-3": ("POR-8-3-090325_S70_L001_R1_001.fastq.gz", "POR-8-3-090325_S70_L001_R2_001.fastq.gz"),
}


def find_genbank_files():
    """Return GenBank file paths on Sockeye."""
    return [
        SOCKEYE_GENBANK / "NC_023133.1_chloroplast.gbk",
        SOCKEYE_GENBANK / "MT483997.1_mitochondrion.gbk",
    ]


def main():
    print("=== Setting up Sockeye agent ===")
    agent_home = SshSource(
        host="sockeye",
        path=Path("/scratch/st-shallam-1/pwy_group/metasmith"),
    ).AsSource()
    smith = Agent(
        home=agent_home,
        runtime=ContainerRuntime.APPTAINER,
        setup_commands=[
            'module load gcc/9.4.0',
            'module load apptainer/1.3.1',
        ],
    )
    print("Agent ready")

    print("\n=== Loading resources & transforms ===")
    base_res = [DataInstanceLibrary.Load(MLIB / "resources/containers")]
    t_transforms = [
        TransformInstanceLibrary.Load(MLIB / "transforms/transcriptomics"),
    ]

    print("\n=== Creating input library ===")
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    inputs_dir = tmp / "inputs.xgdb"
    inputs = DataInstanceLibrary(inputs_dir)
    for tl in ["sequences.yml", "transcriptomics.yml"]:
        inputs.AddTypeLibrary(MLIB / "data_types" / tl)

    # Experiment
    experiment = inputs.AddValue(
        "porphyridium_experiment.txt",
        "porphyridium_transcriptomics",
        "transcriptomics::experiment",
    )

    # GenBank files for organellar genomes (on Sockeye)
    gbk_files = find_genbank_files()
    for gbk_path in gbk_files:
        print(f"  GenBank: {gbk_path.name}")
        inputs.AddItem(gbk_path, "sequences::gbk")

    # Chromosome assembly (on Sockeye)
    inputs.AddItem(SOCKEYE_ASSEMBLY, "sequences::assembly")

    # RNA-seq read pairs
    for sample_name, (r1_file, r2_file) in SAMPLES.items():
        r1_path = SOCKEYE_RAW / r1_file
        r2_path = SOCKEYE_RAW / r2_file
        pair = inputs.AddValue(
            f"{sample_name}_pair.txt", sample_name,
            "sequences::read_pair", parents={experiment},
        )
        inputs.AddItem(r1_path, "sequences::zipped_forward_short_reads", parents={pair})
        inputs.AddItem(r2_path, "sequences::zipped_reverse_short_reads", parents={pair})

    inputs.Save()

    print("\n=== Generating workflow ===")
    targets = TargetBuilder()
    targets.Add("transcriptomics::organellar_gene_count_table")
    targets.Add("transcriptomics::deseq2_results")
    targets.Add("transcriptomics::contamination_report")
    task = smith.GenerateWorkflow(
        samples=list(inputs.AsSamples("transcriptomics::experiment")),
        resources=base_res + [inputs],
        transforms=t_transforms,
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

    # Render DAG
    print("\n=== Rendering workflow DAG ===")
    dag_dir = MLIB / "reports"
    dag_dir.mkdir(parents=True, exist_ok=True)
    try:
        task.plan.RenderDAG(dag_dir / "organellar_workflow_dag.svg")
        task.plan.RenderDAG(dag_dir / "organellar_workflow_dag.png")
        print(f"DAG saved to {dag_dir}")
    except Exception as e:
        print(f"DAG rendering failed (non-fatal): {e}")

    print("\n=== Staging workflow ===")
    smith.StageWorkflow(task, on_exist="update", verify_external_paths=False)
    print("Staged")

    print("\n=== Running workflow (SLURM) ===")
    with open(MLIB / "secrets/slurm_account_sockeye") as f:
        SLURM_ACCOUNT = f.readline().strip()

    smith.RunWorkflow(
        task,
        config_file=smith.GetNxfConfigPresets()["slurm"],
        params=dict(
            slurmAccount=SLURM_ACCOUNT,
            process=dict(scratch="false"),
        ),
    )
    print("Workflow submitted to SLURM")

    print("\n=== Waiting for completion ===")
    results_path = smith.GetResultSource(task).GetPath()
    t0 = time.time()
    timeout = 43200  # 12h
    last_print = 0
    while not (results_path / "_metadata").exists():
        elapsed = time.time() - t0
        if elapsed > timeout:
            print(f"TIMEOUT after {timeout}s")
            sys.exit(1)
        if elapsed - last_print >= 300:
            print(f"  waiting... {elapsed/60:.0f}min elapsed")
            last_print = elapsed
        time.sleep(10)

    smith.CheckWorkflow(task)
    results = DataInstanceLibrary.Load(results_path)
    print(f"\nWorkflow completed in {(time.time()-t0)/60:.0f}min")

    print("\n=== Results ===")
    for path, type_name, endpoint in results.Iterate():
        full_path = path if path.is_absolute() else results_path / path
        print(f"  {type_name}: {full_path} (exists={full_path.exists()})")

    # Cache organellar BAMs to arc storage
    print("\n=== Caching organellar BAMs ===")
    import subprocess
    bam_cache = SOCKEYE_BAM_CACHE
    subprocess.run(
        ["ssh", "sockeye", f"mkdir -p {bam_cache}"],
        check=True,
    )
    for path, type_name, endpoint in results.Iterate():
        if type_name == "transcriptomics::organellar_bam":
            full_path = path if path.is_absolute() else results_path / path
            dest = bam_cache / full_path.name
            print(f"  Caching {full_path.name} -> {dest}")
            subprocess.run(
                ["ssh", "sockeye", f"cp {full_path} {dest}"],
                check=True,
            )

    print("\nDone!")


if __name__ == "__main__":
    main()
