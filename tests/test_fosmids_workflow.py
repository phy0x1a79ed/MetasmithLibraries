"""End-to-end fosmid pipeline test.

Builds a workflow for the fosmid scaffolding + pool-size estimation
chain using the staged fixtures under
``data/metasmith-libraries/staged/fosmids_test/`` and (optionally)
runs it through the configured container runtime.

By default this exercises only workflow generation (matches the
``test_transform_compatibility.py`` style, no containers required).
Set ``FOSMIDS_E2E=1`` to drive the full container execution and
verify outputs end-to-end.
"""
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from metasmith.python_api import (  # noqa: E402
    Agent,
    Source,
    ContainerRuntime,
    DataTypeLibrary,
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
    Resources,
    Size,
    Duration,
)

DOMAINS = [
    "assembly",
    "fosmids",
    "functionalAnnotation",
    "logistics",
    "metabolicModelling",
    "metagenomics",
    "pangenome",
    "responseSurface",
    "transcriptomics",
]

FIXTURES = (
    Path(__file__).resolve().parent.parent.parent.parent.parent
    / "data"
    / "metasmith-libraries"
    / "staged"
    / "fosmids_test"
)


def _require_fixtures():
    expected = [
        "reads.fq.gz",
        "read_metadata.json",
        "host.fna",
        "backbone.fna",
        "fwd_ends.fna",
        "rev_ends.fna",
        "ends_table.json",
        "insert.fna",
    ]
    missing = [n for n in expected if not (FIXTURES / n).exists()]
    if missing:
        pytest.skip(
            f"fosmids fixtures missing: {missing}. "
            f"Run: mamba run -n datascience python {FIXTURES}/generate_fixtures.py"
        )


def _build_inputs(staging: Path) -> tuple[DataInstanceLibrary, DataInstanceLibrary]:
    """Materialize per-sample and per-run libraries for the fosmid pipeline.

    Returns (samples_lib, fosmid_resources_lib).

    - ``samples_lib`` carries the per-sample read inputs. ``read_metadata``
      is the sample-key; ``short_reads_pe`` is parented to it so the
      planner ties them together for each sample.
    - ``fosmid_resources_lib`` carries per-run inputs that aren't tied to a
      single sample: the host genome, vector backbone, and end-sequence
      bundle (forward/reverse fastas parented to the manifest JSON).
    """
    samples_dir = staging / "samples"
    if samples_dir.exists():
        shutil.rmtree(samples_dir)
    samples_lib = DataInstanceLibrary(samples_dir)
    samples_lib.Purge()
    samples_lib.AddTypeLibrary(
        namespace="sequences",
        lib=DataTypeLibrary.Load(ROOT / "data_types/sequences.yml"),
    )
    meta_p = samples_lib.AddItem(
        str(FIXTURES / "read_metadata.json"), "sequences::read_metadata"
    )
    samples_lib.AddItem(
        str(FIXTURES / "reads.fq.gz"),
        "sequences::short_reads_pe",
        parents={meta_p},
    )
    samples_lib.Save()

    res_dir = staging / "fosmid_resources"
    if res_dir.exists():
        shutil.rmtree(res_dir)
    fosmid_res = DataInstanceLibrary(res_dir)
    fosmid_res.Purge()
    fosmid_res.AddTypeLibrary(
        namespace="sequences",
        lib=DataTypeLibrary.Load(ROOT / "data_types/sequences.yml"),
    )
    fosmid_res.AddTypeLibrary(
        namespace="fosmids",
        lib=DataTypeLibrary.Load(ROOT / "data_types/fosmids.yml"),
    )
    fosmid_res.AddItem(str(FIXTURES / "host.fna"), "sequences::background_genome")
    fosmid_res.AddItem(str(FIXTURES / "backbone.fna"), "fosmids::vector_backbone")
    table_p = fosmid_res.AddItem(
        str(FIXTURES / "ends_table.json"), "fosmids::end_sequences_table"
    )
    fosmid_res.AddItem(
        str(FIXTURES / "fwd_ends.fna"),
        "fosmids::end_sequences_forward",
        parents={table_p},
    )
    fosmid_res.AddItem(
        str(FIXTURES / "rev_ends.fna"),
        "fosmids::end_sequences_reverse",
        parents={table_p},
    )
    fosmid_res.Save()
    return samples_lib, fosmid_res


def _make_agent(staging: Path) -> Agent:
    home = Source.FromLocal(staging / "agent_home")
    return Agent(home=home, runtime=ContainerRuntime.DOCKER)


def _wait_for_run(staging: Path, task_key: str, timeout_s: int = 1800) -> Path:
    """Poll the metasmith run log until the workflow completes.

    The launcher backgrounds nextflow via nohup, so RunWorkflow returns
    while the run is still in flight. The run is finished when the
    string ``run completed at`` appears in the latest ``main.log``.
    """
    workspace = staging / "agent_home" / "runs" / task_key
    internals = workspace / "_metasmith"
    deadline = time.time() + timeout_s
    last_log: Path | None = None
    while time.time() < deadline:
        if internals.exists():
            log_dirs = sorted(
                p for p in internals.glob("logs.*") if "latest" not in p.name
            )
            if log_dirs:
                last_log = log_dirs[-1] / "main.log"
                if last_log.exists():
                    text = last_log.read_text(errors="ignore")
                    if "run completed at" in text:
                        return last_log
                    if "ERROR" in text and "nextflow" in text.lower():
                        # nextflow itself reports ERROR lines on failure;
                        # surface the log so pytest shows it
                        raise RuntimeError(
                            f"workflow {task_key} failed; tail of {last_log}:\n"
                            f"{''.join(text.splitlines(keepends=True)[-40:])}"
                        )
        time.sleep(5)
    raise TimeoutError(
        f"workflow {task_key} did not finish within {timeout_s}s "
        f"(latest log: {last_log})"
    )


def _load_resources_and_transforms():
    resources = [
        DataInstanceLibrary.Load(ROOT / f"resources/{n}") for n in ["containers", "lib"]
    ]
    transforms = [
        TransformInstanceLibrary.Load(ROOT / f"transforms/{d}") for d in DOMAINS
    ]
    return resources, transforms


def test_fosmids_scaffold_workflow_generation():
    """Plan a workflow that targets fosmids::scaffolds from raw paired reads.

    Verifies the planner can resolve the full chain: bbduk QC → background
    filter → megahit → fosmids::scaffold. No container execution.
    """
    _require_fixtures()
    staging = Path(tempfile.mkdtemp(prefix="fosmids_gen_"))
    try:
        samples_lib, fosmid_res = _build_inputs(staging)
        resources, transforms = _load_resources_and_transforms()
        targets = TargetBuilder()
        targets.Add("fosmids::scaffolds")

        agent = _make_agent(staging)
        task = agent.GenerateWorkflow(
            samples=samples_lib.AsSamples("sequences::read_metadata"),
            resources=resources + [fosmid_res],
            transforms=transforms,
            targets=targets,
        )
        assert task.ok, "workflow generation failed for fosmids::scaffolds"
        assert len(task.plan.steps) > 0, "no steps in plan"
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def test_fosmids_pool_size_workflow_generation():
    """Plan a workflow that targets fosmids::pool_size_estimate."""
    _require_fixtures()
    staging = Path(tempfile.mkdtemp(prefix="fosmids_pool_"))
    try:
        samples_lib, fosmid_res = _build_inputs(staging)
        resources, transforms = _load_resources_and_transforms()
        targets = TargetBuilder()
        targets.Add("fosmids::pool_size_estimate")

        agent = _make_agent(staging)
        task = agent.GenerateWorkflow(
            samples=samples_lib.AsSamples("sequences::read_metadata"),
            resources=resources + [fosmid_res],
            transforms=transforms,
            targets=targets,
        )
        assert task.ok, "workflow generation failed for fosmids::pool_size_estimate"
        assert len(task.plan.steps) > 0, "no steps in plan"
    finally:
        shutil.rmtree(staging, ignore_errors=True)


@pytest.mark.skipif(
    os.environ.get("FOSMIDS_E2E", "") not in {"1", "true", "yes"},
    reason="full container execution opt-in: set FOSMIDS_E2E=1 to run",
)
def test_fosmids_pipeline_end_to_end():
    """Stage and execute the fosmid pipeline through containers.

    Verifies (a) scaffolds fasta is non-empty and contains the expected
    insert id, (b) pool_size_estimate JSON has size > 0, (c) no stage
    output is empty.
    """
    _require_fixtures()
    if shutil.which("docker") is None:
        pytest.skip("docker runtime not available on PATH")

    staging = Path(tempfile.mkdtemp(prefix="fosmids_run_"))
    try:
        samples_lib, fosmid_res = _build_inputs(staging)
        resources, transforms = _load_resources_and_transforms()

        # target both products in a single workflow run
        targets = TargetBuilder()
        targets.Add("fosmids::scaffolds")
        targets.Add("fosmids::pool_size_estimate")

        agent = _make_agent(staging)
        agent.Deploy()
        task = agent.GenerateWorkflow(
            samples=samples_lib.AsSamples("sequences::read_metadata"),
            resources=resources + [fosmid_res],
            transforms=transforms,
            targets=targets,
        )
        assert task.ok, "workflow generation failed"
        agent.StageWorkflow(task, on_exist="clear")
        # The fixture is tiny; cap every step at minimal resources so a run
        # with all 5 steps in flight (each in its own container) fits on
        # small CI/dev hosts without thrashing.
        agent.RunWorkflow(
            task,
            resource_overrides={
                "all": Resources(
                    cpus=1,
                    memory=Size.GB(1),
                    duration=Duration(hours=2),
                ),
            },
        )
        _wait_for_run(staging, task._key)
        agent.CheckWorkflow(task)

        # walk the run's results directory for produced outputs.
        # metasmith buckets outputs by `<namespace>-<type>/` so we match by
        # parent directory rather than filename.
        results_dir = staging / "agent_home" / "runs" / task._key / "results"
        assert results_dir.exists(), f"results dir not found: {results_dir}"

        def _find_output(parent_name: str, suffix: str) -> Path | None:
            bucket = results_dir / parent_name
            if not bucket.exists():
                return None
            for p in bucket.iterdir():
                if p.is_file() and p.suffix == suffix and p.stat().st_size > 0:
                    return p
            return None

        scaffolds_path = _find_output("fosmids-scaffolds", ".fna")
        pool_path = _find_output("fosmids-pool_size_estimate", ".json")

        assert scaffolds_path is not None, (
            f"scaffolds fasta missing/empty under {results_dir}; "
            f"present: {[p.name for p in results_dir.iterdir()]}"
        )
        assert pool_path is not None, (
            f"pool_size_estimate missing/empty under {results_dir}; "
            f"present: {[p.name for p in results_dir.iterdir()]}"
        )

        # scaffolds fasta should reference each requested insert id
        with open(scaffolds_path) as f:
            scaffold_text = f.read()
        with open(FIXTURES / "ends_table.json") as f:
            ends = json.load(f)
        for iid in ends["insert_ids"]:
            assert iid in scaffold_text, (
                f"insert id [{iid}] missing from scaffolds fasta {scaffolds_path}"
            )

        with open(pool_path) as f:
            pool = json.load(f)
        assert pool.get("size", 0) > 0, f"pool_size <= 0: {pool}"
        # keep workspace on success too — diagnostic for failed runs is via FOSMIDS_E2E_KEEP=1
    finally:
        if os.environ.get("FOSMIDS_E2E_KEEP", "") not in {"1", "true", "yes"}:
            shutil.rmtree(staging, ignore_errors=True)
        else:
            print(f"[FOSMIDS_E2E_KEEP] preserving staging dir: {staging}")
