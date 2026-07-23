#!/usr/bin/env python3
"""Plan the full metagenomics workflow starting from short reads.

Stops at DAG generation. The planner resolves:

  reads --> seqkit_reads --> read_qc_stats
  reads + read_qc_stats --> bbduk --> clean_short_reads
  clean_short_reads --> megahit --> assembly
  assembly --> prodigal --> orfs --> diamond_uniref50 + kofamscan
  assembly --> metabuli                                      (per-contig taxonomy)
  reads + assembly --> assembly_stats --> bam + per-contig/per-bp coverage
       --> metabat2 / semibin2 / comebin --> bin_fasta + contig_to_bin_table
            --> checkm2 --> aggregator --> quality_bin_fasta
                                          --> skani_dedup --> cluster_table
       (one binner's bin_fasta) --> gtdbtk                   (per-bin taxonomy)
  reads --> phyloFlash                                       (SSU rRNA taxonomy)

External DBs (UniRef50, KOFAM, metabuli, GTDB, phyloFlash) are auto-resolved
by including the transforms/logistics library in the plan.

Configuration — set via environment (or edit the defaults):
  MSM_READS_R1/R2   paired short reads to plan from   (REQUIRED)
  MSM_SRC           metasmith source checkout         (optional; else installed)
"""
import os
import sys
from pathlib import Path

# metasmith must be importable; set MSM_SRC to a source checkout if not installed.
if os.environ.get("MSM_SRC"):
    sys.path.insert(0, os.environ["MSM_SRC"])
from metasmith.python_api import (
    Agent, Source, ContainerRuntime,
    DataInstanceLibrary, TransformInstanceLibrary,
    TargetBuilder,
)

R1      = Path(os.environ.get("MSM_READS_R1", "<reads-R1.fq.gz>"))  # paired reads R1
R2      = Path(os.environ.get("MSM_READS_R2", "<reads-R2.fq.gz>"))  # paired reads R2
OUT_DIR = Path("results/metag_workflow")

MLIB = Path(__file__).resolve().parent.parent

inputs = DataInstanceLibrary(OUT_DIR.resolve() / "inputs.xgdb")
for tl in ["sequences.yml", "alignment.yml", "ref.yml", "annotation.yml", "taxonomy.yml", "binning.yml", "binning_local.yml"]:
    inputs.AddTypeLibrary(MLIB / "data_types" / tl)

meta = inputs.AddValue("reads_metadata.json", {"parity": "paired", "length_class": "short"}, "sequences::read_metadata")
pair = inputs.AddValue("read_pair.txt", "sample_1", "sequences::read_pair", parents={meta})
inputs.AddItem(R1.resolve(), "sequences::zipped_forward_short_reads", parents={pair})
inputs.AddItem(R2.resolve(), "sequences::zipped_reverse_short_reads", parents={pair})
inputs.Save()

smith = Agent(home=Source.FromLocal(OUT_DIR.resolve() / "msm_home"), runtime=ContainerRuntime.DOCKER)

targets = TargetBuilder()
targets.Add("sequences::read_qc_stats")
targets.Add("sequences::orfs")
targets.Add("sequences::assembly_stats")
targets.Add("sequences::assembly_per_contig_coverage")
targets.Add("sequences::assembly_per_bp_coverage")
targets.Add("annotation::diamond_uniref50_results")
targets.Add("annotation::kofamscan_results")
targets.Add("taxonomy::metabuli")
targets.Add("taxonomy::phyloflash_summary")
targets.Add("binning_local::cluster_table")

# Per-binner fan-out: distinct TargetSpec parents force a separate
# checkm + gtdbtk instance for each binner's bins (without parent
# constraints the planner would pick one binner to satisfy each).
mb_bin = targets.Add("sequences::metabat2_bin_fasta")
sb_bin = targets.Add("sequences::semibin2_bin_fasta")
cb_bin = targets.Add("sequences::comebin_bin_fasta")
for parent in (mb_bin, sb_bin, cb_bin):
    targets.Add("taxonomy::checkm_stats", parents={parent})
    targets.Add("taxonomy::gtdbtk",       parents={parent})
targets.Add("binning::metabat2_contig_to_bin_table")
targets.Add("binning::semibin2_contig_to_bin_table")
targets.Add("binning::comebin_contig_to_bin_table")

task = smith.GenerateWorkflow(
    samples=list(inputs.AsSamples("sequences::read_metadata")),
    resources=[DataInstanceLibrary.Load(MLIB / "resources" / "env"), inputs],
    transforms=[
        TransformInstanceLibrary.Load(MLIB / "transforms" / "logistics"),
        TransformInstanceLibrary.Load(MLIB / "transforms" / "assembly"),
        TransformInstanceLibrary.Load(MLIB / "transforms" / "metagenomics"),
        TransformInstanceLibrary.Load(MLIB / "transforms" / "functionalAnnotation"),
    ],
    targets=targets,
)

task.plan.RenderDAG(str(OUT_DIR.resolve() / "dag"), format="svg")
