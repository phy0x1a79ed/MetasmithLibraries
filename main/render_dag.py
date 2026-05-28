#!/usr/bin/env python3
"""Render a DAG of all deep-learning transforms using metasmith.

Builds a workflow that targets every output in TARGETS (forces the full
dependency chain incl. shardFasta + saprot's esmfold+foldseek_3di legs),
then calls task.plan.RenderDAG(...) to write an SVG locally.
"""
import sys
from pathlib import Path

import launch_dl_embeddings as L
from metasmith.python_api import DataInstanceLibrary, TransformInstanceLibrary, TargetBuilder

ORFS = L.FOSMIDS_FAA  # smaller input -> same DAG, faster planning

needed_weights = []
seen = set()
for t in L.TARGETS:
    _, _, w_key, _ = L.TARGETS[t]
    if w_key and w_key not in seen:
        needed_weights.append(w_key)
        seen.add(w_key)

inputs = L.build_input_library(ORFS, needed_weights, "dag_all_targets", force_rebuild=True)

containers = DataInstanceLibrary.Load(L.DL_LIB / "resources" / "containers")
transforms = [
    TransformInstanceLibrary.Load(L.DL_LIB / "transforms" / "functionalAnnotation"),
    TransformInstanceLibrary.Load(L.DL_LIB / "transforms" / "logistics"),
]

tb = TargetBuilder()
for t in L.TARGETS:
    type_str, *_ = L.TARGETS[t]
    tb.Add(type_str)

smith = L.get_agent()
print("generating workflow with all targets...")
task = smith.GenerateWorkflow(
    samples=list(inputs.AsSamples("sequences::orfs")),
    resources=[containers, inputs],
    transforms=transforms,
    targets=tb,
)
if not task.ok:
    print(f"ERROR: planning failed: {task.plan}")
    sys.exit(1)

print(f"plan: {len(task.plan.steps)} steps, key={task.GetKey()}")

out_base = Path(__file__).parent / "cache" / "dl_transforms_dag"
out_base.parent.mkdir(parents=True, exist_ok=True)
rendered = task.plan.RenderDAG(str(out_base), format="svg")
svg_path = Path(str(out_base) + ".svg").resolve()
print(f"\nDAG SVG: {svg_path}")
