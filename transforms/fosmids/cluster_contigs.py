"""Dedup the per-pool assembler contigs into a representative ("clustered") set.

Ports scadc `resolve_inserts/cluster.ipynb`: pool every assembly's long contigs
(>= min_contig_len), run an all-vs-all blastn, build a contiguous-pident
(nident/qlen) similarity matrix, and agglomeratively cluster (complete linkage,
distance_threshold = 1 - identity_threshold).  Clustering is done **across all
assemblers together** so near-identical spades/megahit twins of the same clone
collapse into one representative (cross-assembler dedup); `cluster_membership.csv`
records which input contig (and its assembler/pool) fell into each cluster.

Aggregation: a single `fosmids::recovery_experiment` node groups the whole run;
every assembly is parented to it, so `group_by=exp` makes ONE clustering job see
all pools' assemblies via `InputGroup`.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
asm       = model.AddRequirement(lib.GetType("sequences::assembly"), parents={exp})
img_blast = model.AddRequirement(lib.GetType("containers::blast.oci"))
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
out_fa    = model.AddProduct(lib.GetType("fosmids::clustered_contigs"))
out_mem   = model.AddProduct(lib.GetType("fosmids::cluster_membership"))

def protocol(context: ExecutionContext):
    asm_paths = context.InputGroup(asm)
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
    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {prep_script} {asm_arg}",
    )

    # ---------------------------------------------------------------
    # Stage 2: all-vs-all blastn (scadc settings: evalue 1000, perc_identity 50)
    context.ExecWithContainer(
        image=img_blast,
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
    context.ExecWithContainer(
        image=img_pyds,
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
