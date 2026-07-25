"""The ECSPr measurement: fosmid GPR evidence -> conductance effects + significance.

    gpr_table + conditions + atom_pairs + direction_ratios + metabolite_names
        -> ecspr::results

One transform, deliberately. The deployed method splits this across four
(`evidence_weights`, `addition_weights`, `effects`, `significance`) plus a
separate graph-building lane, but every one of those boundaries is an internal
staging decision, not a fork in the method: the intermediate artifacts have a
single producer and a single consumer, and none of them is a thing anyone asks
for on its own. Collapsing them makes the contract state the claim -- this
evidence, under these conditions, on this reference basis, measures this.

The three reference parquets are NETWORK-AGNOSTIC (static functions of the
MNXR/MNXM id space) and so are unpinned shared inputs. `metabolite_names` is not
optional decoration: `atom_pairs` carries metabolite ids and canonical atom
RANKS, so a node of the atom-transfer graph is a (metabolite, rank) key with no
chemical meaning until it is decoded, and the conditions' source/sink hubs
cannot be located on the graph without it.

`conditions` IS per-experiment, and pinned to the experiment for that reason:
the condition set is the measurement's own claim about what it is testing. The
engine has no default set and must not be given one.

LINEAGE. `gpr` is pinned to `exp`, which is the constraint this stage was
waiting on. It forces the GPR table to be the one built from THIS run's
recovered inserts rather than any table the planner could otherwise reach, and
in doing so pulls the whole recovery chain -- reads, host filter, assembly,
junction split, dedup, ORFs, the four annotation lanes -- into the plan behind
it. Drop the pin and the measurement silently becomes a measurement of
something else.

THE NULL IS A SEPARATE PIPELINE and is not an input here. Significance is scored
against draws over a metagenomic ORF pool, size-matched to each unit's ORF
count; that basis is built once, frozen, and shared, and it takes no part in
this experiment's lineage. It arrives as a staged reference when this transform
is implemented -- see `ecsprGround/null.py` in the deployed tree
(`group_by=null_draw_spec`, no `recovery_experiment` requirement).

CONTRACT ONLY -- the protocol is a stub pending the port. The logic to port is
the deployed atom lane: `ecspr/{evidence,addition}_weights.py` and
`ecsprGround/{effects,significance}.py`, with the engine primitives already
staged here as `lib::ecspr_{build,graph,directed}.py`.
"""
from metasmith.python_api import *

lib        = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model      = Transform()

exp        = model.AddRequirement(lib.GetType("fabfos::experiment"))
# THE seam. Pinned to the experiment -- see LINEAGE above.
gpr        = model.AddRequirement(lib.GetType("annotation::gpr_table"), parents={exp})
# The experiment's own claim about what it is testing.
conditions = model.AddRequirement(lib.GetType("ecspr::conditions"), parents={exp})
# The frozen reference basis: shared, unpinned, static in the MNXR/MNXM id space.
pairs      = model.AddRequirement(lib.GetType("ecspr::atom_pairs"))
direction  = model.AddRequirement(lib.GetType("ecspr::direction_ratios"))
names      = model.AddRequirement(lib.GetType("ecspr::metabolite_names"))
img_ecspr  = model.AddRequirement(lib.GetType("env::ecspr.env"))
build      = model.AddRequirement(lib.GetType("lib::ecspr_build.py"))
graphlib   = model.AddRequirement(lib.GetType("lib::ecspr_graph.py"))
directed   = model.AddRequirement(lib.GetType("lib::ecspr_directed.py"))
out        = model.AddProduct(lib.GetType("ecspr::results"))


def protocol(context: ExecutionContext):
    # Intended shape (the four deployed steps, inlined):
    #   1. belief-conservation evidence weights off the GPR table
    #      (per (source, orf, mnxr) E_full/E_dlec/n_orf)
    #   2. per-unit addition maps: contig -> mnxr -> E
    #   3. decode atom_pairs through metabolite_names, apply direction_ratios,
    #      build the base atom-transfer graph, then for each unit solve the
    #      augmented graph and take the signed difference on every condition.
    #      The engine takes NO perturbation argument -- the augmented graph is
    #      rebuilt per unit rather than patched.
    #   4. score each effect against the staged null, size-matched on the unit's
    #      ORF count, and emit the significance table
    raise NotImplementedError(
        "ecspr_measure: contract declared, implementation pending the port of "
        "the deployed atom lane (ecspr/*, ecsprGround/*)"
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(
        cpus=16,
        memory=Size.GB(64),
        duration=Duration(hours=24),
    )
)
