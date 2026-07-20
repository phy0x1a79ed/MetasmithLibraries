"""ECSPr X/Y benchmark SOLVE: frozen benchmark version -> per-shard effect tables.

Thin wrapper around `resources/lib/ecspr_benchmark.py solve`. For every
(facet, element) shard it rebuilds the perturbed graph for each condition in the
answer key and emits the signed conductance change `g_perturbed - g_base` on
every panel edge. The engine's contract -- a difference and not a log ratio,
silent conditions emitted as ZERO rather than dropped -- is documented at length
in that module and is not restated here.

WHY THIS IS A SEPARATE DOMAIN (ecsprBenchmark)
----------------------------------------------
Same reason ecsprNetA and ecsprDirected are separate domains: a loaded domain is
structurally selectable by the planner, so co-loading two producers of one type
makes lane choice a tiebreak instead of an intent. `library.domains_for()` keeps
exactly one of {ecsprNetA, ecsprNetB, ecsprBenchmark} in the domain list, and
that selection -- `--network benchmark` -- is the only thing that puts this lane
in a plan.

WHY NO ANNOTATION LANE RUNS UPSTREAM OF THIS
--------------------------------------------
The benchmark is a FROZEN contract. `ecspr::benchmark_inputs` (the X tree)
already carries the per-reaction evidence weight E for every build, and
`ecspr::benchmark_base_graphs` already carries the induced per-(facet, element)
graphs. So the whole evidence chain -- orfcall, the four annotation lanes,
compile_evidence, evidence_weights, addition_weights, base_graphs -- is upstream
of the freeze, not upstream of this transform. That is the point of the
benchmark: X is fixed, so a score measures the SOLVER, not the annotation.

WHY THE UNIVERSE IS REQUIRED AND NOT OPTIONAL
---------------------------------------------
`--universe` is optional to the CLI but not to us. X deliberately ships no atom
mapping, so a gain-of-function condition -- which must look up the reactions it
inserts and their atom-transit edge weights -- has nowhere to read them from
except `ecspr::benchmark_universe`. Omitting it does not fail fast; it fails on
the first GOF condition, deep into a 16-shard run. Requiring it makes the
planner refuse a tree that lacks it, before anything is scheduled.

WHY ecspr_solver.py IS REQUIRED THOUGH THE PROTOCOL NEVER READS IT
-------------------------------------------------------------------
`ecspr_benchmark.py` does `import ecspr_solver` at module top (after inserting
its own directory on sys.path), so `solve` cannot import without it staged
beside it in the container. Requiring the lib resource stages it into the shared
container dir and the sibling import then resolves. Exact same shape as
ecsprDirected requiring ecspr_directed.py, and as ecsprUndirected requiring
ecspr_solver.py, neither of which context.Input() it either.

THE 16 SHARDS RUN AS ONE TASK, AND THAT IS A KNOWN LIMITATION -- READ THIS
--------------------------------------------------------------------------
The work is 4 facets (netA/netB x iML1515/iECDH10B) x 4 elements (C, N, S, P) =
16 independent shards, and the engine's CLI is deliberately built per-shard
(`--facet`, `--element`, one `--out` each). The natural metasmith expression is
16 sibling task instances. This transform does NOT do that -- it loops over the
shards inside a single task -- and the reason is structural, not laziness:

  * A metasmith transform fans out over the INSTANCES of its requirements. Every
    requirement here (answer key, base graphs, universe, observations, X) is ONE
    staged item covering all 16 shards, so there is exactly one way to satisfy
    this transform and the planner instantiates it exactly once. Nothing in the
    requirement set distinguishes shard from shard.
  * Making the fan-out real needs a per-shard ENDPOINT -- a typed item per
    (facet, element) that the shard id can be read off -- so that 16 distinct
    satisfying assignments exist. That is a change to the benchmark tree's type
    contract (a new `ecspr::benchmark_shard` staged 16 times), not a change to
    this file.
  * And even then it is not sufficient on its own. This project has already been
    bitten by exactly that: an UNCONSTRAINED downstream target makes the planner
    pick ONE producer instead of fanning out over all of them. The fan-out only
    materialises when the downstream consumer (here `merge_benchmark.py`) states
    per-instance `parents={...}` requirements over the shard endpoints. So the
    fix is a matched pair of edits -- new shard type + parents= on merge -- and
    shipping only half of it yields a plan that looks right and silently scores
    one sixteenth of the benchmark.

Rather than emit a plan that quietly solves one shard, the loop is explicit and
in-task: it solves all 16 and the product is the full shard directory. The cost
is scheduling granularity (one long task instead of 16 parallel ones), not
correctness or completeness. Facets are DISCOVERED from the base-graph directory
rather than hardcoded, so a benchmark version that adds or renames a build is
picked up without editing this file; the element list is the engine's own fixed
`choices=ELEMENTS`.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
# NOTE (flagged, not silently decided): unlike ecsprNetA/gem_base_graphs.py this
# transform does NOT require `fosmids::recovery_experiment`. Nothing it consumes
# is per-experiment -- every input is a frozen, version-pinned item of the
# benchmark tree, and the benchmark exists precisely to be independent of any
# one recovery experiment. The consequence is the one that comment describes in
# reverse: the products carry no experiment in their ancestry, so a downstream
# consumer that asks for them with `parents={exp}` will not reach this producer.
# `merge_benchmark.py` accordingly asks without a parents= constraint. If the
# benchmark is ever pulled into an experiment-scoped plan alongside the fosmid
# chain, adding the experiment requirement here is the fix -- and it must be
# added to merge's requirement too, or the pair stops matching.
key       = model.AddRequirement(lib.GetType("ecspr::benchmark_answer_key"))
bases     = model.AddRequirement(lib.GetType("ecspr::benchmark_base_graphs"))
universe  = model.AddRequirement(lib.GetType("ecspr::benchmark_universe"))
obs       = model.AddRequirement(lib.GetType("ecspr::benchmark_observations"))
xinputs   = model.AddRequirement(lib.GetType("ecspr::benchmark_inputs"))
env       = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
bench     = model.AddRequirement(lib.GetType("lib::ecspr_benchmark.py"))
solver    = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
shards    = model.AddProduct(lib.GetType("ecspr::benchmark_observations_solved"))

ELEMENTS = ["C", "N", "S", "P"]

def protocol(context: ExecutionContext):
    ikey  = context.Input(key)
    ibase = context.Input(bases)
    iuni  = context.Input(universe)
    iobs  = context.Input(obs)
    ix    = context.Input(xinputs)
    ibench = context.Input(bench)
    oshards = context.Output(shards)
    # --workers is passed EXPLICITLY rather than left to the engine's
    # $SLURM_CPUS_PER_TASK default. The default is read inside the container,
    # and whether the scheduler's environment survives into the container image
    # is a runtime detail we do not control; `cpus` from the resources line is
    # the same number and is always present. Falling back to 1 keeps the serial
    # reference -- the engine guarantees a parallel run is byte-identical to it,
    # so the fallback changes runtime and nothing else.
    workers = context.params.get("cpus") or 1
    elements = " ".join(ELEMENTS)
    context.ExecWithContainer(
        image=env,
        cmd=f"""for facet_dir in {ibase.container}/*/; do
            facet=$(basename "$facet_dir")
            for el in {elements}; do
                python {ibench.container} solve \
                    --key-dir {ikey.container} \
                    --base-dir {ibase.container} \
                    --universe {iuni.container} \
                    --observations {iobs.container} \
                    --reactions {ix.container}/reactions.tsv \
                    --facet "$facet" --element "$el" \
                    --workers {workers} \
                    --out {oshards.container}/"$facet"/"$el".tsv
            done
        done""",
    )
    return ExecutionResult(
        manifest=[{shards: oshards.local}],
        success=oshards.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=key,
    # cpus=32: the batch size this lane is provisioned for. It is NOT decorative
    # and it is NOT a BLAS thread count -- the engine pins OMP/OPENBLAS/MKL/
    # NUMEXPR to 1 before importing numpy (oversubscription measured >10x slower
    # in this project's null builder) and parallelises with a FORK POOL over
    # CONDITIONS instead. So 32 here means 32 worker processes each running whole
    # single-threaded solves, which is the shape that actually scales. It is also
    # what `--workers` receives above.
    resources=Resources(cpus=32, memory=Size.GB(64), duration=Duration(hours=12)),
)
