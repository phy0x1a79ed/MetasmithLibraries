#!/usr/bin/env python3
"""Plan a PathoLogic (Pathway Tools) workflow for a nucleotide assembly FASTA.

Targets both ptools outputs (cyc folder + parsed CSV tables); the planner
fans out orfs -> {deepec, kofamscan, diamond_uniref50} -> ptools_annotation_gather
-> pathologic. Stops at DAG generation.

Configuration — set via environment (or edit the defaults):
  MSM_ASSEMBLY   nucleotide assembly FASTA to annotate   (REQUIRED)
  MSM_SRC        metasmith source checkout to import      (optional; else use
                                                           an installed metasmith)
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

ASSEMBLY = Path(os.environ.get("MSM_ASSEMBLY", "<assembly.fna>"))  # nucleotide assembly FASTA
OUT_DIR = Path("results/pathologic")

MLIB = Path(__file__).resolve().parent.parent

inputs = DataInstanceLibrary(OUT_DIR.resolve() / "inputs.xgdb")
inputs.AddTypeLibrary(MLIB / "data_types" / "sequences.yml")
inputs.AddItem(ASSEMBLY.resolve(), "sequences::assembly")
inputs.Save()

smith = Agent(home=Source.FromLocal(OUT_DIR.resolve() / "msm_home"), runtime=ContainerRuntime.DOCKER)

targets = TargetBuilder()
targets.Add("annotation::pgdb_archive")
targets.Add("annotation::pgdb_csv_tables")

task = smith.GenerateWorkflow(
    samples=list(inputs.AsSamples("sequences::assembly")),
    resources=[DataInstanceLibrary.Load(MLIB / "resources" / "env"), inputs],
    transforms=[
        TransformInstanceLibrary.Load(MLIB / "transforms" / "logistics"),
        TransformInstanceLibrary.Load(MLIB / "transforms" / "metagenomics"),
        TransformInstanceLibrary.Load(MLIB / "transforms" / "functionalAnnotation"),
    ],
    targets=targets,
)

task.plan.RenderDAG(str(OUT_DIR.resolve() / "dag"), format="svg")
