"""ECSPr X/Y benchmark MERGE: per-shard effect tables -> the one observations.tsv.

Thin wrapper around `resources/lib/ecspr_benchmark.py merge`. It globs the shard
tree, concatenates every non-empty table, sorts on (facet, condition_id,
edge_id) and writes the single table the scorer consumes.

WHY THIS IS A TRANSFORM AT ALL, AND NOT A TAIL OF THE SOLVE
------------------------------------------------------------
Two reasons, both about what a re-run costs.

First, the merge is the join point. The moment the shards become 16 real task
instances (see the long note in `solve_benchmark.py` -- this is deliberately
written so that change touches only requirements, not the protocol), the merge is
the node that has to see all of them. Folding it into the solve would mean the
16-way split could not be made without also rewriting the concatenation.

Second, cache granularity. Merging is pure concatenation over hashed shard
inputs; separating it means re-deriving the merged table is seconds of work and
never re-solves. A solve that ended in a write of its own merged copy would put
a 12-hour task behind every change to how the shards are assembled.

WHY THE ENGINE REFUSES AN EMPTY MERGE, AND WHY WE LET IT
--------------------------------------------------------
`cmd_merge` returns 1 when no shard produced rows, rather than writing an empty
observations.tsv. That refusal is the thing standing between "the benchmark
scored zero conditions" and "the benchmark reports nothing and the scorer reads
a well-formed empty table". The container command is not wrapped in any `|| true`
for exactly that reason: a non-zero exit must fail the task.

Note that an EMPTY SHARD is not the same as an empty merge and is not an error --
`cmd_solve` writes a header-only table for a (facet, element) with no panel edges
or conditions on that element, and merge skips it. Only ALL shards being empty
trips the refusal.

WHY ecspr_solver.py IS REQUIRED THOUGH MERGE NEVER SOLVES ANYTHING
-------------------------------------------------------------------
`ecspr_benchmark.py` imports `ecspr_solver` at MODULE top, before argparse ever
sees which subcommand was asked for, so `merge` cannot import without it staged
beside it in the container any more than `solve` can. Requiring the lib resource
stages it into the shared container dir; the sibling import then resolves.
"""
from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
# NO parents= constraint on the shards, and that is a decision worth naming.
# `solve_benchmark.py` requires no experiment, so its product carries no
# experiment in its ancestry and a `parents={exp}` requirement here would not
# match it -- the planner would report this target unsatisfiable rather than
# name the missing edge (the failure mode ecsprNetA/gem_base_graphs.py documents
# from the other side).
#
# THE FAN-OUT CAVEAT, restated here because this is the file it lands in: today
# `shards` is ONE directory holding all 16 (facet, element) tables, so one
# unconstrained requirement is correct and complete. If the solve is ever split
# into 16 task instances, this single requirement becomes WRONG in a way that
# raises nothing -- an unconstrained target lets the planner satisfy it with ONE
# producer instead of all sixteen, and the merge then silently concatenates one
# shard into a table that looks like a finished benchmark. The split therefore
# requires per-instance `parents={...}` requirements HERE, added in the same
# change; it is not a solve-side edit.
sharddir = model.AddRequirement(lib.GetType("ecspr::benchmark_observations_solved"))
env      = model.AddRequirement(lib.GetType("containers::ecspr.oci"))
bench    = model.AddRequirement(lib.GetType("lib::ecspr_benchmark.py"))
solver   = model.AddRequirement(lib.GetType("lib::ecspr_solver.py"))
result   = model.AddProduct(lib.GetType("ecspr::benchmark_result"))

def protocol(context: ExecutionContext):
    ishards = context.Input(sharddir)
    ibench  = context.Input(bench)
    oresult = context.Output(result)
    context.ExecWithContainer(
        image=env,
        cmd=f"""python {ibench.container} merge \
            --shard-dir {ishards.container} \
            --out {oresult.container}""",
    )
    return ExecutionResult(
        manifest=[{result: oresult.local}],
        success=oresult.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=sharddir,
    # cpus=1. This is pandas.concat plus a sort over a table of tens of
    # thousands of rows -- there is no thread flag on `merge`, nothing to widen,
    # and the engine pins the maths libraries to one thread regardless. Asking
    # for more would reserve cores and idle them behind a scheduler queue for a
    # step measured in seconds.
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=1)),
)
