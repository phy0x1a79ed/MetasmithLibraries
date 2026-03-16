#!/usr/bin/env python
"""Generate the Metasmith metabolomics workflow DAG using render_dag().

Uses the actual transforms (jgi_loader, differential_analysis, pathway_enrichment)
and their declared data types to produce an accurate pipeline DAG.

Run: conda run -n msm_env python3 results/generate_dag.py
"""
import sys
import tempfile
from pathlib import Path

MLIB = Path(__file__).resolve().parent.parent
sys.path.insert(0, "/home/tony/workspace/msm/Metasmith/src")

from metasmith.python_api import (
    Agent, ContainerRuntime, Source,
    DataInstanceLibrary, TransformInstanceLibrary,
    TargetBuilder,
)

# Local agent (planning only, no SSH needed)
tmp = Path(tempfile.mkdtemp())
agent_home = Source.FromLocal(tmp)
smith = Agent(home=agent_home, runtime=ContainerRuntime.APPTAINER)

# Transform library — metabolicModelling transforms
transforms = [
    TransformInstanceLibrary.Load(MLIB / "transforms/metabolicModelling"),
]

# Container library — metabolomics-python type is defined in transform-local types
containers = DataInstanceLibrary(tmp / "containers.xgdb")
containers.AddTypeLibrary(MLIB / "data_types" / "containers.yml")
# Also load transform-local type defs (metabolomics-python.oci is defined there)
containers.AddTypeLibrary(
    MLIB / "transforms/metabolicModelling/_metadata/types/containers.yml"
)
containers.AddItem(
    MLIB / "resources/containers/metabolomics-python.oci",
    "containers::metabolomics-python.oci",
)
containers.Save()

# Input library — mock the JGI dataset input
inputs_dir = tmp / "inputs.xgdb"
inputs = DataInstanceLibrary(inputs_dir)
for tl in ["metabolomics.yml", "containers.yml"]:
    inputs.AddTypeLibrary(MLIB / "data_types" / tl)


def mock(name, is_dir=False):
    p = tmp / name
    if is_dir:
        p.mkdir(exist_ok=True)
    else:
        p.touch()
    return p


# JGI metabolomics dataset (starting input)
inputs.AddItem(mock("jgi_data", is_dir=True), "metabolomics::jgi_metabolomics_dataset")

inputs.Save()

# Targets: pathway enrichment (terminal output of the 3-step pipeline)
targets = TargetBuilder()
targets.Add("metabolomics::metabolomics_pathway_enrichment")

print("Generating workflow plan...", flush=True)
task = smith.GenerateWorkflow(
    samples=list(inputs.AsSamples("metabolomics::jgi_metabolomics_dataset")),
    resources=[containers, inputs],
    transforms=transforms,
    targets=targets,
)

if not task.ok:
    print(f"FAILED to generate plan: {task.plan}")
    sys.exit(1)

steps = task.plan.steps
print(f"\nPlan OK — {len(steps)} steps\n")
for step in steps:
    name = Path(step.transform._path).stem
    prods = [i.dtype_name for g in step.produces for i in g]
    print(f"  Step {step.order}: {name} -> {prods}")

# Render DAG — ensure graphviz dot is on PATH
import os
_env_bin = Path(sys.executable).parent
os.environ["PATH"] = f"{_env_bin}:{os.environ.get('PATH', '')}"

out_dir = MLIB / "results" / "report_package01"
out_dir.mkdir(parents=True, exist_ok=True)

out_svg = out_dir / "workflow_dag.svg"
task.plan.RenderDAG(out_svg, blacklist_namespaces={"lib", "containers"})
print(f"\nDAG written to: {out_svg}")

out_png = out_dir / "workflow_dag.png"
task.plan.RenderDAG(out_png, blacklist_namespaces={"lib", "containers"})
print(f"DAG written to: {out_png}")
