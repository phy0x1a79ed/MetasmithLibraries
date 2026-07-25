"""Map pCC1 vector-backbone junctions in one pool's contigs and cut at them.

One step, two products. Blast the fosmid vector backbone against every contig
this pool's assemblers produced, record where the backbone lands, and emit the
pieces cut at those positions. Two kinds of hit matter and the map distinguishes
them:

* **terminal** -- backbone at (or within a tolerance of) a contig end. That is a
  real clone boundary: the assembly ran off the insert into the vector.
* **internal** -- backbone mid-contig. Two clones the assembler joined through
  their shared backbone, i.e. a chimera, and the position is where to cut.

Mapping and cutting are deliberately ONE transform rather than two. The split is
a pure function of the map plus the contigs the map was called against, so a
separate step would only re-stage the same two assemblies to re-derive an
in-memory interval list -- and would need the same host-filtered pin on the same
assemblies to guarantee it cut the contigs the junctions were called on. The
`junction_map` CSV is still emitted as its own product, so the evidence is
inspectable and the split stays auditable; it is simply not a staging boundary.

This replaces the coverage-matrix chimera call. Detection is now positive
evidence (the vector is *there*) rather than an inference from a depth
discontinuity, so it does not need the per-experiment coverage matrix and runs
per pool, in that pool's own assembly coordinates. There is also deliberately
**no trimming step** any more: the old chain split on a coverage discontinuity
and then trimmed ragged low-depth ends, both inferences from the depth profile.
A junction is positive evidence of a boundary, so the cut is the whole
operation.

Every output carries the header-baked interval provenance the domain already
uses (`source=<contig_id> start= end= action=`), so a downstream reader can
always recover which input contig a piece came from and what was done to it.

Aggregation: `group_by=meta` -- one job per pool, seeing that pool's megahit and
spades contigs together, so a junction called in one assembler's copy of a clone
is recorded alongside the other's.

Host-filter coercion lives here. This is the first consumer of the assemblies in
the recovery chain (it used to be `cluster_contigs`), so pinning the consumed
assembly's lineage to `host_filtered_short_reads` is what forces
`reads -> background_filter -> {megahit, spades}` upstream. Move it and host
filtering silently drops out of the plan.

CONTRACT ONLY -- the protocol is a stub pending the implementation migration.
The chimera-calling and interval/header logic this supersedes is in git at
`chimera_split.py` / `coverage_trim.py` (both removed).
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
# The experiment node is the root of the run. Both per-run references below hang
# off it, so a plan cannot silently reach for some other run's backbone or host.
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
meta      = model.AddRequirement(lib.GetType("sequences::read_metadata"), parents={exp})
# Lineage constraint only -- the protocol never reads it. See the note above:
# pinning the assemblies to a host-filtered ancestor is what pulls
# `background_filter` into the plan, since `host_filtered_short_reads` extends
# `clean_short_reads` and the assemblers accept it unchanged.
hf        = model.AddRequirement(lib.GetType("sequences::host_filtered_short_reads"), parents={meta})
# Both assemblers explicitly: a generic `sequences::assembly` requirement lets
# the planner satisfy the DAG with one of them, but the intended pipeline runs
# megahit AND spades per pool and carries both through to the dedup.
asm_mh    = model.AddRequirement(lib.GetType("sequences::megahit_assembly"), parents={meta, hf})
asm_sp    = model.AddRequirement(lib.GetType("sequences::spades_assembly"), parents={meta, hf})
backbone  = model.AddRequirement(lib.GetType("fosmids::vector_backbone"), parents={exp})
img_blast = model.AddRequirement(lib.GetType("env::blast.env"))
img_pyds  = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
out_jm    = model.AddProduct(lib.GetType("fosmids::junction_map"))
out_split = model.AddProduct(lib.GetType("fosmids::split_contigs"))


def protocol(context: ExecutionContext):
    # Intended shape:
    #   1. makeblastdb over this pool's contigs (both assemblers, union)
    #   2. blastn the backbone against it
    #   3. classify each hit terminal|internal against a contig-end tolerance,
    #      emit one row per hit with the derived junction_pos -> out_jm
    #   4. for each contig, cut at every internal junction_pos and drop the
    #      backbone span itself at terminal junctions; write each piece with
    #      source=/start=/end=/action= in the header -> out_split
    raise NotImplementedError(
        "junction_split: contract declared, implementation pending the migration pass"
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=meta,
    resources=Resources(
        cpus=4,
        memory=Size.GB(16),
        duration=Duration(hours=2),
    )
)
