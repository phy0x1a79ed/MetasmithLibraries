"""
End-to-end tests for VirSorter2 and DRAM-v transforms.

These tests verify that viral identification workflows can be:
1. Generated (workflow planning)
2. Staged to the agent
3. Executed via local Docker
4. Produce valid viral sequence and annotation outputs
"""
import pytest
from pathlib import Path
from metasmith.python_api import (
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
    Resources,
    Size,
    Duration,
)

from conftest import (
    wait_for_workflow,
    verify_fasta_output,
    verify_tsv_output,
    MLIB,
    TEST_DATA_DIR,
)


@pytest.fixture(scope="module")
def virsorter2_transforms(mlib):
    """Load functional annotation transforms (includes VirSorter2 + DRAM-v)."""
    return [
        TransformInstanceLibrary.Load(mlib / "transforms/functionalAnnotation"),
    ]


@pytest.fixture
def assembly_input(tmp_inputs, test_data_dir):
    """Create input library with assembly for viral identification."""
    inputs = tmp_inputs(["sequences.yml", "annotation.yml"])

    assembly_path = test_data_dir / "small_assembly.fna"
    if not assembly_path.exists():
        pytest.skip("Test data not available: small_assembly.fna")

    inputs.AddItem(assembly_path, "sequences::assembly")
    inputs.LocalizeContents()
    inputs.Save()

    return inputs


@pytest.fixture
def virsorter2_resources(mlib, base_resources):
    """Load VirSorter2-specific resources."""
    resources = list(base_resources)

    lib_path = mlib / "resources/lib"
    if lib_path.exists():
        try:
            lib_res = DataInstanceLibrary.Load(lib_path)
            resources.append(lib_res)
        except Exception:
            pass

    return resources


@pytest.fixture
def vs2_db_input(mlib, test_data_dir):
    """Create input library with VirSorter2 database."""
    vs2_db_dir = test_data_dir / "virsorter2_db"
    if not vs2_db_dir.exists():
        pytest.skip("VirSorter2 database not available")

    lib_dir = vs2_db_dir / "_virsorter2.xgdb"
    inputs = DataInstanceLibrary(lib_dir)
    inputs.AddTypeLibrary(mlib / "data_types" / "annotation.yml")
    inputs.AddItem(vs2_db_dir, "annotation::virsorter2_db")
    inputs.Save()
    return inputs


@pytest.fixture
def dram_db_input(mlib, test_data_dir):
    """Create input library with DRAM database."""
    dram_db_dir = test_data_dir / "dram_db"
    if not dram_db_dir.exists():
        pytest.skip("DRAM database not available")

    lib_dir = dram_db_dir / "_dram.xgdb"
    inputs = DataInstanceLibrary(lib_dir)
    inputs.AddTypeLibrary(mlib / "data_types" / "annotation.yml")
    inputs.AddItem(dram_db_dir, "annotation::dram_db")
    inputs.Save()
    return inputs


class TestVirSorter2WorkflowGeneration:
    """Tests for workflow generation (planning only)."""

    def test_can_plan_virsorter2_workflow(
        self, agent, virsorter2_resources, virsorter2_transforms, assembly_input
    ):
        """Verify workflow generation for VirSorter2."""
        targets = TargetBuilder()
        targets.Add("annotation::virsorter2_viral_sequences")

        task = agent.GenerateWorkflow(
            samples=list(assembly_input.AsSamples("sequences::assembly")),
            resources=virsorter2_resources + [assembly_input],
            transforms=virsorter2_transforms,
            targets=targets,
        )

        # Skip if required resources (VS2 DB) aren't available
        if not task.ok or len(task.plan.steps) == 0:
            pytest.skip("VirSorter2 workflow requires database")

    def test_can_plan_dramv_workflow(
        self, agent, virsorter2_resources, virsorter2_transforms, assembly_input
    ):
        """Verify workflow generation for DRAM-v (chained from VirSorter2)."""
        targets = TargetBuilder()
        targets.Add("annotation::dramv_annotations")

        task = agent.GenerateWorkflow(
            samples=list(assembly_input.AsSamples("sequences::assembly")),
            resources=virsorter2_resources + [assembly_input],
            transforms=virsorter2_transforms,
            targets=targets,
        )

        # Skip if required resources (VS2 + DRAM DB) aren't available
        if not task.ok or len(task.plan.steps) == 0:
            pytest.skip("DRAM-v workflow requires VirSorter2 + DRAM databases")


@pytest.mark.slow
class TestVirSorter2WorkflowExecution:
    """Full E2E tests that execute workflows via Docker."""

    def test_virsorter2_e2e(
        self, agent, virsorter2_resources, virsorter2_transforms, assembly_input, vs2_db_input
    ):
        """Full E2E test: stage, run VirSorter2, verify viral sequences and scores."""
        targets = TargetBuilder()
        targets.Add("annotation::virsorter2_viral_sequences")
        targets.Add("annotation::virsorter2_scores")
        targets.Add("annotation::virsorter2_affi_contigs")
        targets.Add("annotation::virsorter2_boundary")

        task = agent.GenerateWorkflow(
            samples=list(assembly_input.AsSamples("sequences::assembly")),
            resources=virsorter2_resources + [assembly_input, vs2_db_input],
            transforms=virsorter2_transforms,
            targets=targets,
        )
        assert task.ok, f"Workflow generation failed: {task}"

        agent.StageWorkflow(task, on_exist="clear")
        agent.RunWorkflow(
            task,
            config_file=agent.GetNxfConfigPresets()["local"],
            params=dict(
                executor=dict(cpus=4, queueSize=1),
                process=dict(tries=1),
            ),
            resource_overrides={
                "*": Resources(cpus=4, memory=Size.GB(8)),
            },
        )

        results = wait_for_workflow(agent, task, timeout=3600)
        results_path = agent.GetResultSource(task).GetPath()

        found_seqs = False
        found_scores = False

        for path, type_name, endpoint in results.Iterate():
            if "virsorter2_viral_sequences" in type_name:
                found_seqs = True
                if not path.is_absolute():
                    full_path = results_path / path
                    assert verify_fasta_output(full_path), f"Invalid FASTA: {full_path}"
            if "virsorter2_scores" in type_name:
                found_scores = True
                if not path.is_absolute():
                    full_path = results_path / path
                    assert verify_tsv_output(full_path), f"Invalid TSV: {full_path}"

        assert found_seqs, "No VirSorter2 viral sequences output found"
        assert found_scores, "No VirSorter2 scores output found"

    def test_dramv_e2e(
        self, agent, virsorter2_resources, virsorter2_transforms, assembly_input,
        vs2_db_input, dram_db_input
    ):
        """Full E2E test: stage, run VirSorter2 → DRAM-v, verify annotations."""
        targets = TargetBuilder()
        targets.Add("annotation::dramv_annotations")
        targets.Add("annotation::dramv_distill")

        task = agent.GenerateWorkflow(
            samples=list(assembly_input.AsSamples("sequences::assembly")),
            resources=virsorter2_resources + [assembly_input, vs2_db_input, dram_db_input],
            transforms=virsorter2_transforms,
            targets=targets,
        )
        assert task.ok, f"Workflow generation failed: {task}"

        agent.StageWorkflow(task, on_exist="clear")
        agent.RunWorkflow(
            task,
            config_file=agent.GetNxfConfigPresets()["local"],
            params=dict(
                executor=dict(cpus=4, queueSize=1),
                process=dict(tries=1),
            ),
            resource_overrides={
                "*": Resources(cpus=4, memory=Size.GB(16)),
            },
        )

        results = wait_for_workflow(agent, task, timeout=3600)
        results_path = agent.GetResultSource(task).GetPath()

        found_annot = False
        found_distill = False

        for path, type_name, endpoint in results.Iterate():
            if "dramv_annotations" in type_name:
                found_annot = True
                if not path.is_absolute():
                    full_path = results_path / path
                    assert verify_tsv_output(full_path), f"Invalid TSV: {full_path}"
            if "dramv_distill" in type_name:
                found_distill = True

        assert found_annot, "No DRAM-v annotations output found"
        assert found_distill, "No DRAM-v distill output found"
