# MetasmithLibraries

Hallam Lab transform, data-type, and container library for the
[Metasmith](https://github.com/hallamlab/metasmith) workflow engine.

Metasmith plans and executes data workflows by composing typed transforms.
This repo is the registry of types, container images, and transform
implementations that planner draws from.

## Contents

```
data_types/      # YAML type definitions (15 namespaces)
resources/
  containers/    # 62 OCI container references (.oci → docker:// URI)
  lib/           # Reusable analysis scripts (clustering, plotting, etc.)
transforms/      # Transform implementations grouped by domain
  amplicon/              #  2 transforms
  assembly/              # 13 transforms
  functionalAnnotation/  # 27 transforms
  logistics/             # 23 transforms (downloaders, format conversion)
  metabolicModelling/    #  5 transforms
  metagenomics/          #  9 transforms (+ binning/, taxonomy/)
  pangenome/             #  2 transforms
  responseSurface/       #  1 transform
  transcriptomics/       # 19 transforms
tests/           # Pytest workflow & E2E tests
envs/            # Conda env for running tests outside Metasmith
analysis/        # Downstream analysis code (not part of the library)
main/            # Operator scripts (DAG rendering, planner probes, run launchers)
```

Data types (`data_types/*.yml`) cover: `alignment`, `amplicon`, `annotation`,
`binning`, `clustering`, `containers`, `lib`, `media_optimization`,
`metabolomics`, `ncbi`, `pangenome`, `ref`, `sequences`, `taxonomy`,
`transcriptomics`.

Transform domains include genome/metagenome **assembly** (flye, hifiasm-meta,
megahit, metaspades, …), **binning** (metabat2, semibin2, comebin, checkm),
**taxonomy** (gtdbtk, fastani, kraken2/bracken, sylph, …), **functional
annotation** (eggnog-mapper, kofamscan, interproscan, antismash, bakta,
dram, busco, diamond, deepec, …), **protein embeddings & structure**
(esmfold, esm-c, prott5, ankh, saprot, foldseek), **transcriptomics**
(star, salmon, stringtie, braker3, deseq2/pydeseq2), and **logistics**
(weight/DB downloaders, SRA/assembly fetchers, read format conversion,
sharding, container prefetch).

## Building the library

The build compiles every `data_types/*.yml`, every `resources/*/` directory,
and every `transforms/*/` directory into the `_metadata/` indexes that
Metasmith's planner loads at runtime.

```bash
./dev.sh -b
```

The script delegates to `msm build all --types … --uniques … --transforms …`.
It uses `msm` if it's on `PATH`; otherwise it falls back to
`../Metasmith/dev.sh -r build all …`, so you need either:

- a Metasmith install with `msm` on `PATH` (matching the dev-branch CLI signature), or
- a sibling checkout of [hallamlab/Metasmith](https://github.com/hallamlab/metasmith) plus a Python env satisfying `envs/base.yml` (and `envs/dev.yml` for tests).

After the build, `_metadata/` directories under `resources/` and
`transforms/` reflect the current type graph. They're committed to the
repo so consumers can pull the library without rebuilding.

### Adding a container

1. Add the type to `data_types/containers.yml` (`extends: container`,
   declare what it `provides`).
2. Create `resources/containers/<name>.oci` containing the registry URI
   (e.g. `docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0`).
3. Rerun `./dev.sh -b`.

The build fails if any `.oci` file lacks a matching type definition.

### Adding a transform

```python
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
inp   = model.AddRequirement(lib.GetType("sequences::assembly"))
out   = model.AddProduct(lib.GetType("annotation::eggnog_table"))

def protocol(context: ExecutionContext):
    in_path  = context.Input(inp)
    out_path = context.Output(out)
    context.ExecWithContainer(
        image=lib.GetResource("containers::eggnog-mapper.oci"),
        cmd=f"emapper.py -i {in_path.container} -o {out_path.container}",
    )
    return ExecutionResult(manifest=[{out: out_path.local}], success=out_path.local.exists())

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=inp,
    resources=Resources(cpus=8, memory=Size.GB(32), duration=Duration(hours=12)),
)
```

See `transforms/_template.py` for the minimal skeleton.

Every type referenced via `lib.GetType("namespace::type")` must exist in
`data_types/<namespace>.yml` or the build fails.

## Tests

Tests use the fixtures defined in `tests/conftest.py` (`agent`,
`base_resources`, `mlib`, `tmp_inputs`).

```bash
pytest tests/                                       # full suite
pytest tests/test_binning_workflow.py -v            # one workflow
pytest tests/ -m "not slow"                         # skip E2E
./dev.sh --test-comebin                             # named E2E shortcut
```

Workflow-planning tests call `agent.GenerateWorkflow(...)` and assert on
the resulting plan. End-to-end tests additionally call `StageWorkflow()` /
`RunWorkflow()` and wait for completion; they're marked `@pytest.mark.slow`
and need a configured execution backend.

## Layout conventions

- `transforms/*/_metadata/` is build-generated — do not hand-edit.
- `data_types/`, `resources/*/_metadata/index.yml` (for declaring containers),
  and the transform `.py` files themselves are the hand-edited surface.
- Disabled transforms live in `transforms/*/_disabled/` or are renamed
  `<name>.py.disabled` so the build skips them.

## Related repos

- [hallamlab/metasmith](https://github.com/hallamlab/metasmith) — the planner / executor
- This library is consumed by Metasmith agents at workflow-generation time;
  it has no Python package of its own.
