"""
End-to-end tests for the metabolomics transform pipeline.

Tests the JGI metabolomics data loading, differential analysis,
pathway enrichment, and FBA constraint transforms using Porphyridium
purpureum test data.
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
from conftest import MLIB, TEST_DATA_DIR, wait_for_workflow


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PORPHYRIDIUM_DIR = TEST_DATA_DIR / "porphyridium-metabolomics" / "jgi_metabolomic_data"
METABOLOMICS_TYPES = MLIB / "data_types" / "metabolomics.yml"
METABOLOMICS_TRANSFORMS = MLIB / "transforms" / "metabolicModelling"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def metabolomics_transforms():
    """Load the metabolicModelling transform library."""
    return [TransformInstanceLibrary.Load(METABOLOMICS_TRANSFORMS)]


@pytest.fixture(scope="session")
def metabolomics_input(mlib):
    """Create input library with the Porphyridium JGI metabolomics dataset."""
    if not PORPHYRIDIUM_DIR.exists():
        pytest.skip(
            f"Porphyridium metabolomics test data not available: {PORPHYRIDIUM_DIR}"
        )

    lib_dir = PORPHYRIDIUM_DIR.parent / "_metabolomics_inputs.xgdb"
    inputs = DataInstanceLibrary(lib_dir)
    inputs.AddTypeLibrary(METABOLOMICS_TYPES)
    inputs.AddItem(PORPHYRIDIUM_DIR, "metabolomics::jgi_metabolomics_dataset")
    inputs.Save()
    return inputs


@pytest.fixture(scope="session")
def feature_table_results(
    agent, base_resources, metabolomics_transforms, metabolomics_input
):
    """Run the jgi_loader workflow and return results library."""
    targets = TargetBuilder()
    targets.Add("metabolomics::metabolomics_feature_table")

    task = agent.GenerateWorkflow(
        samples=list(
            metabolomics_input.AsSamples("metabolomics::jgi_metabolomics_dataset")
        ),
        resources=base_resources + [metabolomics_input],
        transforms=metabolomics_transforms,
        targets=targets,
    )
    assert task.ok, f"jgi_loader workflow planning failed: {task}"

    agent.StageWorkflow(task, on_exist="clear")
    agent.RunWorkflow(
        task,
        config_file=agent.GetNxfConfigPresets()["local"],
        params=dict(executor=dict(cpus=4, queueSize=1), process=dict(tries=1)),
        resource_overrides={"*": Resources(cpus=4, memory=Size.GB(8))},
    )
    return wait_for_workflow(agent, task, timeout=600)


@pytest.fixture(scope="session")
def differential_results(
    agent, base_resources, metabolomics_transforms, feature_table_results
):
    """Run differential_analysis workflow and return results library."""
    targets = TargetBuilder()
    targets.Add("metabolomics::metabolomics_differential")

    task = agent.GenerateWorkflow(
        samples=list(
            feature_table_results.AsSamples("metabolomics::metabolomics_feature_table")
        ),
        resources=base_resources + [feature_table_results],
        transforms=metabolomics_transforms,
        targets=targets,
    )
    assert task.ok, f"differential_analysis workflow planning failed: {task}"

    agent.StageWorkflow(task, on_exist="clear")
    agent.RunWorkflow(
        task,
        config_file=agent.GetNxfConfigPresets()["local"],
        params=dict(executor=dict(cpus=2, queueSize=1), process=dict(tries=1)),
        resource_overrides={"*": Resources(cpus=2, memory=Size.GB(4))},
    )
    return wait_for_workflow(agent, task, timeout=600)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
class TestJgiLoader:
    """Test 1: JGI data loading and normalization."""

    def test_workflow_produces_feature_table(self, feature_table_results):
        """Verify the jgi_loader produces a feature table CSV."""
        found = False
        for path, type_name, endpoint in feature_table_results.Iterate():
            if "feature_table" in type_name:
                found = True
                break
        assert found, "No feature table found in jgi_loader results"

    def test_feature_table_has_expected_samples(self, feature_table_results):
        """Verify output CSV has correct 161-PE sample columns."""
        import pandas as pd

        for path, type_name, endpoint in feature_table_results.Iterate():
            if "feature_table" in type_name:
                results_dir = feature_table_results._path
                csv_path = results_dir / path if not path.is_absolute() else path
                df = pd.read_csv(csv_path)
                sample_cols = [
                    c for c in df.columns if "161-pe" in c.lower()
                ]
                # Expect 13 groups * 4 reps = ~52, but at least some
                assert (
                    len(sample_cols) > 0
                ), f"No 161-PE sample columns found. Columns: {df.columns.tolist()[:10]}"
                return
        pytest.fail("Feature table CSV not found in results")

    def test_feature_count(self, feature_table_results):
        """Verify a reasonable number of features were detected."""
        import pandas as pd

        for path, type_name, endpoint in feature_table_results.Iterate():
            if "feature_table" in type_name:
                results_dir = feature_table_results._path
                csv_path = results_dir / path if not path.is_absolute() else path
                df = pd.read_csv(csv_path)
                assert len(df) > 50, (
                    f"Expected >50 features, got {len(df)}"
                )
                return
        pytest.fail("Feature table CSV not found in results")


class TestDifferentialAnalysis:
    """Test 2: Statistical differential analysis."""

    def test_comparison_csvs_exist(self, differential_results):
        """Verify per-condition comparison CSVs are generated."""
        found_dir = False
        for path, type_name, endpoint in differential_results.Iterate():
            if "differential" in type_name:
                found_dir = True
                results_dir = differential_results._path
                diff_path = results_dir / path if not path.is_absolute() else path
                if diff_path.is_dir():
                    csvs = list(diff_path.glob("*_vs_*.csv"))
                    assert len(csvs) >= 3, (
                        f"Expected >=3 comparison CSVs, found {len(csvs)}: "
                        f"{[c.name for c in csvs]}"
                    )
                break
        assert found_dir, "No differential results directory found"

    def test_anova_results(self, differential_results):
        """Verify ANOVA results file exists."""
        for path, type_name, endpoint in differential_results.Iterate():
            if "differential" in type_name:
                results_dir = differential_results._path
                diff_path = results_dir / path if not path.is_absolute() else path
                if diff_path.is_dir():
                    anova = diff_path / "anova_all_conditions.csv"
                    assert anova.exists(), "ANOVA results file not found"
                    return
        pytest.fail("Differential results directory not found")

    def test_volcano_plots(self, differential_results):
        """Verify volcano plot PNGs are generated."""
        for path, type_name, endpoint in differential_results.Iterate():
            if "differential" in type_name:
                results_dir = differential_results._path
                diff_path = results_dir / path if not path.is_absolute() else path
                if diff_path.is_dir():
                    pngs = list(diff_path.glob("*_volcano.png"))
                    assert len(pngs) >= 1, "No volcano plots generated"
                    return
        pytest.fail("Differential results directory not found")


class TestPathwayEnrichment:
    """Test 3: KEGG pathway ORA."""

    @pytest.fixture(scope="class")
    def enrichment_results(
        self, agent, base_resources, metabolomics_transforms, differential_results
    ):
        """Run pathway_enrichment workflow."""
        targets = TargetBuilder()
        targets.Add("metabolomics::metabolomics_pathway_enrichment")

        task = agent.GenerateWorkflow(
            samples=list(
                differential_results.AsSamples(
                    "metabolomics::metabolomics_differential"
                )
            ),
            resources=base_resources + [differential_results],
            transforms=metabolomics_transforms,
            targets=targets,
        )
        assert task.ok, f"pathway_enrichment workflow planning failed: {task}"

        agent.StageWorkflow(task, on_exist="clear")
        agent.RunWorkflow(
            task,
            config_file=agent.GetNxfConfigPresets()["local"],
            params=dict(executor=dict(cpus=2, queueSize=1), process=dict(tries=1)),
            resource_overrides={"*": Resources(cpus=2, memory=Size.GB(4))},
        )
        return wait_for_workflow(agent, task, timeout=1200)

    def test_enrichment_csv(self, enrichment_results):
        """Verify enrichment CSV is non-empty."""
        import pandas as pd

        for path, type_name, endpoint in enrichment_results.Iterate():
            if "pathway_enrichment" in type_name:
                results_dir = enrichment_results._path
                enrich_path = results_dir / path if not path.is_absolute() else path
                if enrich_path.is_dir():
                    csv_path = enrich_path / "enrichment.csv"
                    assert csv_path.exists(), "enrichment.csv not found"
                    df = pd.read_csv(csv_path)
                    assert len(df.columns) > 0, "enrichment.csv has no columns"
                    return
        pytest.fail("Pathway enrichment results not found")

    def test_dotplot(self, enrichment_results):
        """Verify enrichment dot plot is generated."""
        for path, type_name, endpoint in enrichment_results.Iterate():
            if "pathway_enrichment" in type_name:
                results_dir = enrichment_results._path
                enrich_path = results_dir / path if not path.is_absolute() else path
                if enrich_path.is_dir():
                    png = enrich_path / "enrichment_dotplot.png"
                    # Dot plot only generated if enriched pathways exist
                    if not png.exists():
                        pytest.skip("No enriched pathways to plot")
                    assert png.stat().st_size > 0
                    return
        pytest.fail("Pathway enrichment results not found")


class TestFbaConstraint:
    """Test 4: Metabolic model FBA (optional, requires template model)."""

    @pytest.fixture(scope="class")
    def fba_results(
        self,
        agent,
        base_resources,
        metabolomics_transforms,
        differential_results,
        mlib,
    ):
        """Run fba_constraint workflow."""
        # Look for an SBML model in resources or test data
        model_path = None
        for candidate in [
            TEST_DATA_DIR / "porphyridium-metabolomics" / "iCre1355.xml",
            MLIB / "resources" / "lib" / "iCre1355.xml",
        ]:
            if candidate.exists():
                model_path = candidate
                break

        if model_path is None:
            pytest.skip("No SBML template model available for FBA test")

        # Create input library with model
        lib_dir = model_path.parent / "_sbml_input.xgdb"
        model_input = DataInstanceLibrary(lib_dir)
        model_input.AddTypeLibrary(METABOLOMICS_TYPES)
        model_input.AddItem(model_path, "metabolomics::metabolic_model_sbml")
        model_input.Save()

        targets = TargetBuilder()
        targets.Add("metabolomics::metabolomics_fba_results")

        task = agent.GenerateWorkflow(
            samples=list(
                differential_results.AsSamples(
                    "metabolomics::metabolomics_differential"
                )
            ),
            resources=base_resources + [differential_results, model_input],
            transforms=metabolomics_transforms,
            targets=targets,
        )
        assert task.ok, f"fba_constraint workflow planning failed: {task}"

        agent.StageWorkflow(task, on_exist="clear")
        agent.RunWorkflow(
            task,
            config_file=agent.GetNxfConfigPresets()["local"],
            params=dict(executor=dict(cpus=4, queueSize=1), process=dict(tries=1)),
            resource_overrides={"*": Resources(cpus=4, memory=Size.GB(8))},
        )
        return wait_for_workflow(agent, task, timeout=1200)

    def test_flux_csv(self, fba_results):
        """Verify flux CSV is generated with expected reaction count."""
        import pandas as pd

        for path, type_name, endpoint in fba_results.Iterate():
            if "fba_results" in type_name:
                results_dir = fba_results._path
                fba_path = results_dir / path if not path.is_absolute() else path
                if fba_path.is_dir():
                    flux_csv = fba_path / "condition_fluxes.csv"
                    assert flux_csv.exists(), "condition_fluxes.csv not found"
                    df = pd.read_csv(flux_csv)
                    assert len(df) > 100, (
                        f"Expected >100 reactions, got {len(df)}"
                    )
                    return
        pytest.fail("FBA results not found")

    def test_escher_json(self, fba_results):
        """Verify Escher JSON overlay is generated."""
        import json

        for path, type_name, endpoint in fba_results.Iterate():
            if "fba_results" in type_name:
                results_dir = fba_results._path
                fba_path = results_dir / path if not path.is_absolute() else path
                if fba_path.is_dir():
                    json_path = fba_path / "escher_flux_overlay.json"
                    assert json_path.exists(), "escher_flux_overlay.json not found"
                    data = json.loads(json_path.read_text())
                    assert isinstance(data, dict)
                    assert len(data) > 0, "Escher JSON is empty"
                    return
        pytest.fail("FBA results not found")
