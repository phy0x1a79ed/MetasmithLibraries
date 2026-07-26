"""Dedup the junction-split contigs into the final putative insert set.

Ports scadc `resolve_inserts/cluster.ipynb`: pool every input's long contigs
(>= min_contig_len), run an all-vs-all blastn, build a contiguous-pident
(nident/qlen) similarity matrix, and agglomeratively cluster (complete linkage,
distance_threshold = 1 - identity_threshold).  Clustering is done **across all
assemblers together** so near-identical spades/megahit twins of the same clone
collapse into one representative (cross-assembler dedup); `cluster_membership.csv`
records which input contig (and its assembler/pool) fell into each cluster.

This step now runs **last**, not first. Junction calling and splitting moved
ahead of it (see `junction_map.py` / `split_at_junctions.py`), because those work
per pool in that pool's own assembly coordinates. So the input here is the split
pieces rather than raw assemblies, and the representative set that comes out of
the dedup *is* the final insert set -- there is no separate trim afterwards, and
no intermediate `clustered_contigs` type any more.

The cost of splitting first is that the same chimera assembled by both megahit
and spades gets cut twice; that is exactly what this dedup collapses.

Aggregation: a single `fabfos::experiment` node groups the whole run;
each pool's `read_metadata` -- and so its split contigs -- descends from it, so
`group_by=exp` makes ONE dedup job see every pool's pieces via `InputGroup`.

CONTRACT REWIRED, protocol not yet migrated: the body below still reads two
assembly inputs and writes `clustered_contigs`. It needs to read the split
contigs and write the inserts + report instead.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fabfos::experiment"))
# The host-filter coercion that used to live here has moved to `junction_map`,
# which is now the first consumer of the assemblies. Everything reaching this
# step is already downstream of it.
split     = model.AddRequirement(lib.GetType("fabfos::split_contigs"), parents={exp})
img_blast = model.AddRequirement(lib.GetType("env::blast.env"))
img_pyds  = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
out_ins   = model.AddProduct(lib.GetType("fabfos::putative_inserts"))
out_rep   = model.AddProduct(lib.GetType("fabfos::putative_insert_report"))
out_mem   = model.AddProduct(lib.GetType("fabfos::cluster_membership"))

def protocol(context: ExecutionContext):
    # STALE -- kept verbatim for the migration pass. The contract above now takes
    # `split` and produces inserts + report + membership; the body below still
    # reads two assemblies and writes `clustered_contigs`. Guarded so it cannot
    # run half-migrated and quietly produce the wrong thing.
    raise NotImplementedError(
        "cluster_contigs: contract rewired to consume fabfos::split_contigs and "
        "produce putative_inserts + report; body not yet migrated"
    )

    # union of both assemblers' per-pool contigs (cross-assembler dedup input)
    asm_paths = context.InputGroup(asm_mh) + context.InputGroup(asm_sp)
    ofa       = context.Output(out_fa)
    omem      = context.Output(out_mem)

    min_contig_len = context.params.get("min_contig_length", 10000)
    threshold      = context.params.get("cluster_identity", 0.99)
    threads        = context.params.get("cpus")
    threads_arg    = "" if threads is None else f"-num_threads {threads}"

    # ---------------------------------------------------------------
    # Stage 1: pool contigs >= min_contig_len, rekey stably, record (assembler,
    # source_pool) parsed from each file's `<pool>__<asm>` stem.
    asm_arg = " ".join(str(p.container) for p in asm_paths)
    prep_script = "prep_contigs.py"
    with open(prep_script, "w") as f:
        f.write(f"""\
import sys, json
from pathlib import Path
from Bio import SeqIO

MIN_LEN = {min_contig_len}
seqs = {{}}
meta = {{}}
i = 0
for path in sys.argv[1:]:
    stem = Path(path).name
    for suf in (".fna", ".fasta", ".fa"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    if "__" in stem:
        pool, assembler = stem.split("__", 1)
    else:
        pool, assembler = stem, "unknown"
    for rec in SeqIO.parse(path, "fasta"):
        if len(rec.seq) < MIN_LEN:
            continue
        i += 1
        k = f"C{{i:05d}}"
        seqs[k] = str(rec.seq)
        meta[k] = dict(assembler=assembler, source_pool=pool,
                       description=rec.description, length=len(rec.seq))

assert len(seqs) > 0, f"no contigs >= {{MIN_LEN}}bp across the staged assemblies"
with open("pooled.fna", "w") as out:
    for k, s in seqs.items():
        out.write(f">{{k}}\\n{{s}}\\n")
with open("contig_meta.json", "w") as j:
    json.dump(meta, j)
print(f"pooled_contigs={{len(seqs)}}")
""")
    context.ExecWithEnv().ifContainerDo(
        env=img_pyds,
        cmd=f"python {prep_script} {asm_arg}",
    )

    # ---------------------------------------------------------------
    # Stage 2: all-vs-all blastn (scadc settings: evalue 1000, perc_identity 50)
    context.ExecWithEnv().ifContainerDo(
        env=img_blast,
        cmd=f"""\
            makeblastdb -dbtype nucl -in pooled.fna -out pooled_db >makeblastdb.log 2>&1
            blastn -evalue 1000 -perc_identity 50 {threads_arg} \
                -query pooled.fna \
                -db pooled_db \
                -outfmt "6 qseqid sseqid qstart qend nident qlen slen" \
                -out blast.tsv >blastn.log 2>&1
        """,
    )

    # ---------------------------------------------------------------
    # Stage 3: cluster within each assembler group, union representatives.
    cluster_script = "cluster.py"
    with open(cluster_script, "w") as f:
        f.write(f"""\
import sys, json
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering

THRESHOLD = {threshold}
ofa  = sys.argv[1]
omem = sys.argv[2]

with open("contig_meta.json") as j:
    meta = json.load(j)
seqs = {{}}
cur = None
with open("pooled.fna") as fa:
    for line in fa:
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
seqs = {{k: "".join(v) for k, v in seqs.items()}}

# contiguous pident per ordered pair: max(nident/qlen)
COLS = "qseqid sseqid qstart qend nident qlen slen".split()
try:
    dfe = pd.read_csv("blast.tsv", sep="\\t", header=None, names=COLS)
except pd.errors.EmptyDataError:
    dfe = pd.DataFrame(columns=COLS)
pident = {{}}
for _, r in dfe.iterrows():
    k = (r.qseqid, r.sseqid)
    pident[k] = max(pident.get(k, 0.0), r.nident / r.qlen if r.qlen else 0.0)

def cluster_group(labels):
    if len(labels) == 1:
        return {{labels[0]: labels[0]}}   # member -> centroid
    c2i = {{c: i for i, c in enumerate(labels)}}
    n = len(labels)
    sim = np.zeros((n, n))
    for (q, s), v in pident.items():
        if q in c2i and s in c2i:
            i, j = c2i[q], c2i[s]
            sim[i, j] = max(sim[i, j], v)
            sim[j, i] = max(sim[j, i], v)
    np.fill_diagonal(sim, 1.0)  # self-identity (Ns from end-mapping otherwise < 1)
    m = AgglomerativeClustering(metric="precomputed", linkage="complete",
                                n_clusters=None, distance_threshold=1 - THRESHOLD)
    m.fit(1 - sim)
    clusters = {{}}
    for lab, c in zip(m.labels_, labels):
        clusters.setdefault(int(lab), []).append(c)
    member2centroid = {{}}
    for members in clusters.values():
        # centroid = max summed similarity to the rest of its cluster
        best, best_score = members[0], -1.0
        for a in members:
            score = sum(sim[c2i[a], c2i[b]] for b in members)
            if score > best_score:
                best, best_score = a, score
        for mm in members:
            member2centroid[mm] = best
    return member2centroid

# cluster ALL contigs together (across assemblers) so spades/megahit twins merge
member2centroid = cluster_group(list(meta.keys()))

representatives = sorted(set(member2centroid.values()))
with open(ofa, "w") as out:
    for c in representatives:
        mv = meta[c]
        out.write(f">{{c}} assembler={{mv['assembler']}} source_pool={{mv['source_pool']}} "
                  f"length={{mv['length']}}\\n{{seqs[c]}}\\n")

rows = []
for member, centroid in sorted(member2centroid.items()):
    mv = meta[member]
    rows.append(dict(cluster=centroid, centroid=centroid, member=member,
                     assembler=mv["assembler"], source_pool=mv["source_pool"],
                     length=mv["length"]))
pd.DataFrame(rows, columns=["cluster", "centroid", "member", "assembler",
                            "source_pool", "length"]).to_csv(omem, index=False)
print(f"representatives={{len(representatives)}} from {{len(meta)}} contigs "
      f"(cross-assembler clustering)")
""")
    context.ExecWithEnv().ifContainerDo(
        env=img_pyds,
        cmd=f"python {cluster_script} {ofa.container} {omem.container}",
    )

    return ExecutionResult(
        manifest=[{out_fa: ofa.local, out_mem: omem.local}],
        success=ofa.local.exists() and omem.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(
        cpus=8,
        memory=Size.GB(32),
        duration=Duration(hours=12),
    )
)
