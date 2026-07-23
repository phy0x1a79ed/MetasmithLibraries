# MetasmithLibraries — Development Notes

## Project structure

```
data_types/          # YAML type definitions (sequences, env, etc.)
resources/           # Data instance libraries (env declarations, reference DBs)
transforms/          # Transform implementations grouped by domain
  logistics/         # Data retrieval & format conversion
  assembly/          # Genome/metagenome assembly
  metagenomics/      # Binning, taxonomy, etc.
  functionalAnnotation/
  amplicon/
  pangenome/
tests/               # Pytest-based workflow & E2E tests
  test_data/         # Test datasets (ORA files, mock reads, assemblies)
  conftest.py        # Shared fixtures: agent, base_resources, tmp_inputs, etc.
```

## Build system

- **Rebuild metadata:** `./dev.sh -b` (requires `msm` CLI from the `msm_env` conda environment)
- Alternatively: `conda run -n msm_env msm build --types data_types --uniques resources/* --transforms transforms/*`
- The build regenerates all `_metadata/` directories from source YAML + transform Python files
- Every env resource file in `resources/env/` **must** have a matching type definition in `data_types/env.yml`, otherwise the build fails
- Every type referenced via `lib.GetType("namespace::type")` in transforms must exist in the corresponding `data_types/*.yml`

## Environments (containers + conda)

A tool's environment is declared generically. `resources/env/<tool>.env` is a
YAML file with an optional `container:` (a `docker://…` OCI URI, used by the
DOCKER/APPTAINER runtimes) and/or an optional `conda:` (a conda env name, used by
the MAMBA runtime). The engine selects which by the single global runtime (see
`ResolveEnvImage` / `GetContainerModel` in the metasmith engine). One `.env`
serves both a containerized run and a mamba run.

### Adding a new env

1. Add the type to `data_types/env.yml` with `extends: env` and a `provides` list
2. Create `resources/env/<name>.env` with `container: docker://…` and/or `conda: <name>`
3. If it has a `conda:` env, add a recipe `envs/tools/<name>.yml` (or re-run
   `python envs/gen_tool_envs.py`, which derives biocontainers specs automatically)
4. Run `./dev.sh -b` to rebuild metadata; `./dev.sh --create-envs` creates the conda test envs

## Writing transforms

- Transforms are Python files using `from metasmith.python_api import *`
- `TransformInstanceLibrary.ResolveParentLibrary(__file__)` loads types from the parent library's `_metadata/`
- Use `model.AddRequirement()` for inputs, `model.AddProduct()` for outputs
- For paired-end reads, use a grouping parent (e.g. `read_pair`) and set `parents={pair}` on both R1/R2 requirements
- `group_by=` in `TransformInstance()` controls how inputs are matched/grouped
- `context.ExecWithContainer(image=, cmd=)` runs commands in the resolved env (container or `mamba run`)
- Container paths: `context.Input(x).container` (path inside container), `.local` (path on host), `.external` (path from outside container)

## Writing tests

- Tests use `conftest.py` fixtures: `agent`, `base_resources`, `mlib`, `tmp_inputs`
- `tmp_inputs(["sequences.yml", ...])` creates a temporary `DataInstanceLibrary` with specified type libraries
- Use `inputs.AddItem(path, type)` for files, `inputs.AddValue(name, dict, type)` for JSON values
- Workflow generation tests: `agent.GenerateWorkflow(samples=, resources=, transforms=, targets=)`
- E2E tests: additionally call `agent.StageWorkflow()`, `agent.RunWorkflow()`, then `wait_for_workflow()`
- Mark E2E tests with `@pytest.mark.slow`
- Run tests with: `conda run -n msm_env pytest tests/<file>.py -k "<pattern>" -v`

## Driver scripts (`main/`)

Top-level scripts that exercise the library end-to-end (build inputs, plan, optionally stage/run).

- `main/diamond_uniref50_from_assembly.py` — minimal planning-only driver: one input type (`sequences::assembly`), one target (`annotation::diamond_uniref50_results`). The canonical example of the simplest possible DAG-generating shape.
- `main/metag_workflow_from_reads.py` — full metagenomics workflow starting from short reads: seqkit_reads → bbduk → megahit → prodigal → {diamond_uniref50, kofamscan}; metabuli; assembly_stats → 3 binners (metabat2/semibin2/comebin) → checkm2 → aggregator → skani_dedup; gtdbtk; phyloFlash. Targets every intermediate (read_qc_stats, orfs, assembly stats/coverage, per-binner contig_to_bin tables, checkm_stats, …) so they all appear in the rendered DAG.
- `main/{diamond_uniref50_from_assembly,metag_workflow_from_reads}_sockeye.py` — HPC (Sockeye-style: SSH + APPTAINER + allocation-coded `/scratch`) end-to-end ports of the two drivers above. Reference inputs/DBs by REMOTE paths, supply pre-staged DBs as resources to prune the `local`-labeled `download*` steps, and run Deploy → Generate → Stage → Run → Wait. The metag one is a **W0/W1 pair**: run `main/metag_setup_sockeye.py --run` first (prefetch the 17 tool containers via `containers::pulled_container`, upload reads, verify DBs), then `main/metag_workflow_from_reads_sockeye.py --run`.

Conventions:

- **Public-repo-safe config** — no hardcoded absolute paths, allocations, usernames, or DB paths. All site-specific values come from env vars with `<placeholder>` defaults (`MSM_SRC`, `MSM_ASSEMBLY`, `MSM_READS_R1/R2`, `MSM_HPC_HOST`, `MSM_SLURM_ACCOUNT`, `MSM_REF_DB_DIR`, …); `MSM_SRC` optionally prepends a metasmith source checkout to `sys.path` (else an installed metasmith is used). HPC drivers `require_configured()` and exit if a placeholder is unfilled.
- Use `inputs.AddValue(name, dict, "namespace::type")` for inline values (e.g. `read_metadata`) and `inputs.AddItem(path, type, parents={...})` to chain lineage (meta → reads → … is the canonical metag pattern when starting from reads)
- Include `transforms/logistics` so the planner auto-resolves external DBs (UniRef50, KOFAM, metabuli, GTDB, phyloFlash) — or supply them as pre-staged resources (`ref::uniref50_diamond_db`, `ref::kofamscan_profiles`/`_ko_list`, `ref::metabuli_ref`, `ref::gtdb`, `ref::phyloflash_db`) to skip the (login-node-only) downloads on HPC
- To force all sibling transforms (e.g. all 3 binners), target a downstream that requires them all (`binning_local::cluster_table` pulls metabat2 + semibin2 + comebin + per-binner checkm + aggregator + skani_dedup), or list each sibling's product as a target
- End with `task.plan.RenderDAG(out, format="svg")` for a planning-only driver

Longer drivers (`launch_dl_embeddings.py`, `probe_planner.py`, `render_dag.py`) keep the same structure but add CLI parsing, remote SshSource agents, and full stage/run/wait flows.

## Conda environment

- Use `msm_env` for running `msm build` and `pytest`
- `conda run -n msm_env <command>` or `conda activate msm_env`
