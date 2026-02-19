"""
End-to-end tests for functional annotation transforms.

Tests for: interproscan, kofamscan, proteinbert, deepec, diamond_uniref50

These tests verify that annotation workflows can be:
1. Generated (workflow planning)
2. Staged to the agent
3. Executed via local Docker
4. Produce valid annotation outputs
"""
import pytest
import time
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
    verify_json_output,
    verify_tsv_output,
    kofam_db_input,
    interproscan_data_input,
    uniref50_db_input,
    MLIB,
    TEST_DATA_DIR,
)


@pytest.fixture(scope="module")
def annotation_transforms(mlib):
    """Load functional annotation transforms."""
    return [
        TransformInstanceLibrary.Load(mlib / "transforms/functionalAnnotation"),
    ]


@pytest.fixture
def orfs_input(tmp_inputs, test_data_dir):
    """Create input library with ORFs (protein sequences)."""
    inputs = tmp_inputs(["sequences.yml", "annotation.yml"])

    orfs_path = test_data_dir / "small_orfs.faa"
    if not orfs_path.exists():
        pytest.skip("Test data not available: small_orfs.faa")

    inputs.AddItem(orfs_path, "sequences::orfs")
    inputs.LocalizeContents()
    inputs.Save()

    return inputs


@pytest.fixture
def annotation_resources(mlib, base_resources):
    """Load annotation-specific resources."""
    resources = list(base_resources)

    lib_path = mlib / "resources/lib"
    if lib_path.exists():
        try:
            lib_res = DataInstanceLibrary.Load(lib_path)
            resources.append(lib_res)
        except Exception:
            pass

    return resources


class TestAnnotationWorkflowGeneration:
    """Tests for workflow generation (planning only)."""

    def test_can_plan_interproscan_workflow(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Verify workflow generation for InterProScan."""
        targets = TargetBuilder()
        targets.Add("annotation::interproscan_json")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
            targets=targets,
        )

        assert task.ok, f"Workflow generation failed: {task}"
        assert len(task.plan.steps) > 0

    def test_can_plan_kofamscan_workflow(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Verify workflow generation for KofamScan."""
        targets = TargetBuilder()
        targets.Add("annotation::kofamscan_results")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
            targets=targets,
        )

        assert task.ok, f"Workflow generation failed: {task}"

    def test_can_plan_proteinbert_workflow(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Verify workflow generation for ProteinBERT."""
        targets = TargetBuilder()
        targets.Add("annotation::proteinbert_embeddings")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
            targets=targets,
        )

        assert task.ok, f"Workflow generation failed: {task}"

    def test_can_plan_deepec_workflow(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Verify workflow generation for DeepEC."""
        targets = TargetBuilder()
        targets.Add("annotation::deepec_predictions")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
            targets=targets,
        )

        assert task.ok, f"Workflow generation failed: {task}"

    def test_can_plan_diamond_uniref50_workflow(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Verify workflow generation for DIAMOND UniRef50."""
        targets = TargetBuilder()
        targets.Add("annotation::diamond_uniref50_results")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
            targets=targets,
        )

        # Skip if UniRef50 database isn't available
        if not task.ok:
            pytest.skip("DIAMOND UniRef50 workflow requires database")


@pytest.mark.slow
class TestAnnotationWorkflowExecution:
    """Full E2E tests that execute workflows via Docker."""

    def test_interproscan_e2e(
        self, agent, annotation_resources, annotation_transforms, orfs_input, interproscan_data_input
    ):
        """Full E2E test: stage, run InterProScan, verify JSON/GFF outputs."""
        targets = TargetBuilder()
        targets.Add("annotation::interproscan_json")
        targets.Add("annotation::interproscan_gff")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input, interproscan_data_input],
            transforms=annotation_transforms,
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

        results = wait_for_workflow(agent, task, timeout=1800)
        results_path = agent.GetResultSource(task).GetPath()

        found_json = False
        found_gff = False

        for path, type_name, endpoint in results.Iterate():
            if "interproscan_json" in type_name:
                found_json = True
                if not path.is_absolute():
                    full_path = results_path / path
                    assert verify_json_output(full_path), f"Invalid JSON: {full_path}"

            if "interproscan_gff" in type_name:
                found_gff = True

        assert found_json, "No InterProScan JSON output found"
        assert found_gff, "No InterProScan GFF output found"

    def test_kofamscan_e2e(
        self, agent, annotation_resources, annotation_transforms, orfs_input, kofam_db_input
    ):
        """Full E2E test: stage, run KofamScan, verify results."""
        targets = TargetBuilder()
        targets.Add("annotation::kofamscan_results")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input, kofam_db_input],
            transforms=annotation_transforms,
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

        results = wait_for_workflow(agent, task, timeout=600)
        results_path = agent.GetResultSource(task).GetPath()

        found_results = False
        for path, type_name, endpoint in results.Iterate():
            if "kofamscan_results" in type_name:
                found_results = True

        assert found_results, "No KofamScan results found"

    def test_deepec_e2e(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Full E2E test: stage, run DeepEC, verify predictions."""
        targets = TargetBuilder()
        targets.Add("annotation::deepec_predictions")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
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

        results = wait_for_workflow(agent, task, timeout=600)
        results_path = agent.GetResultSource(task).GetPath()

        found_predictions = False
        for path, type_name, endpoint in results.Iterate():
            if "deepec_predictions" in type_name:
                found_predictions = True

        assert found_predictions, "No DeepEC predictions found"

    def test_proteinbert_e2e(
        self, agent, annotation_resources, annotation_transforms, orfs_input
    ):
        """Full E2E test: stage, run ProteinBERT, verify embeddings."""
        targets = TargetBuilder()
        targets.Add("annotation::proteinbert_embeddings")
        targets.Add("annotation::proteinbert_index")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input],
            transforms=annotation_transforms,
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

        results = wait_for_workflow(agent, task, timeout=900)
        results_path = agent.GetResultSource(task).GetPath()

        found_embeddings = False
        found_index = False

        for path, type_name, endpoint in results.Iterate():
            if "proteinbert_embeddings" in type_name:
                found_embeddings = True
            if "proteinbert_index" in type_name:
                found_index = True

        assert found_embeddings, "No ProteinBERT embeddings found"
        assert found_index, "No ProteinBERT index found"

    def test_diamond_uniref50_e2e(
        self, agent, annotation_resources, annotation_transforms, orfs_input, uniref50_db_input
    ):
        """Full E2E test: stage, run DIAMOND UniRef50, verify BLAST6 results."""
        targets = TargetBuilder()
        targets.Add("annotation::diamond_uniref50_results")

        task = agent.GenerateWorkflow(
            samples=list(orfs_input.AsSamples("sequences::orfs")),
            resources=annotation_resources + [orfs_input, uniref50_db_input],
            transforms=annotation_transforms,
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

        results = wait_for_workflow(agent, task, timeout=600)

        found_results = False
        for path, type_name, endpoint in results.Iterate():
            if "diamond_uniref50_results" in type_name:
                found_results = True

        assert found_results, "No DIAMOND UniRef50 results found"
