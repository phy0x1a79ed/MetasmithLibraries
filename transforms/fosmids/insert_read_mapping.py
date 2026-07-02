"""Map the pooled reads back to the reference inserts -> BAM + coverage table.

The analysis-side counterpart to the recovery-side `pool_coverage`: instead of
mapping to the clustered contigs to drive splitting/trimming, this maps the pooled
reads to the FINAL reference inserts (the 132 >= 29 kb fosmids) to produce the
"read mapping data" deliverable.

For each pool's reads (staged one file per pool, `<pool>.fq`) it runs
minimap2 -a -x sr --secondary=no -> samtools sort, then:
  * merges all per-pool sorted BAMs into one `alignment::bam` (+ `alignment::csi`)
    -- a single coordinate-sorted alignment of every pool against the inserts;
  * folds each pool's bedgraph into a per-(insert, pool) coverage summary
    (mean/median depth, breadth at 1x/5x) plus an aggregate "all"-pool row,
    emitted as `fosmids::insert_coverage`.

Aggregation mirrors pool_coverage: one `fosmids::recovery_experiment` node groups the
run; the reference inserts and every pool's reads descend from it, so group_by=exp
gives ONE job the single reference fasta plus all pools' reads via InputGroup.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
reference = model.AddRequirement(lib.GetType("fosmids::reference_inserts"), parents={exp})
reads     = model.AddRequirement(lib.GetType("sequences::clean_short_reads"), parents={exp})
img_mm2   = model.AddRequirement(lib.GetType("containers::minimap2.oci"))
img_sam   = model.AddRequirement(lib.GetType("containers::samtools.oci"))
img_bed   = model.AddRequirement(lib.GetType("containers::bedtools.oci"))
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
out_bam   = model.AddProduct(lib.GetType("alignment::bam"))
out_csi   = model.AddProduct(lib.GetType("alignment::csi"))
out_cov   = model.AddProduct(lib.GetType("fosmids::insert_coverage"))

def protocol(context: ExecutionContext):
    iref        = context.Input(reference)
    reads_paths = context.InputGroup(reads)
    obam        = context.Output(out_bam)
    ocsi        = context.Output(out_csi)
    ocov        = context.Output(out_cov)

    threads     = context.params.get("cpus")
    threads_mm2 = "" if threads is None else f"-t {threads}"
    threads_sam = "" if threads is None else f"-@ {threads}"

    pool_bams = []
    for rp in reads_paths:
        pool = rp.local.name
        for suf in ("_unmapped", ".fq", ".fastq", ".gz", ".fasta", ".fna"):
            pool = pool.replace(suf, "")
        context.ExecWithContainer(
            image=img_mm2,
            cmd=f"""\
                minimap2 -a -x sr {threads_mm2} --secondary=no \
                    {iref.container} {rp.container} > aln.{pool}.sam
            """,
        )
        context.ExecWithContainer(
            image=img_sam,
            cmd=f"""\
                samtools sort {threads_sam} -o aln.{pool}.bam aln.{pool}.sam
                rm -f aln.{pool}.sam
            """,
        )
        context.ExecWithContainer(
            image=img_bed,
            cmd=f"""\
                bedtools genomecov -ibam aln.{pool}.bam -bga > cov__{pool}.bedgraph
            """,
        )
        pool_bams.append(f"aln.{pool}.bam")

    # merge all per-pool sorted BAMs into one coordinate-sorted alignment + CSI index
    context.ExecWithContainer(
        image=img_sam,
        cmd=f"""\
            samtools merge {threads_sam} -f {obam.container} {' '.join(pool_bams)}
            samtools index -c {threads_sam} {obam.container} {ocsi.container}
            rm -f {' '.join(pool_bams)}
        """,
    )

    # fold each pool's bedgraph into a per-(insert, pool) coverage summary
    cov_script = "insert_coverage.py"
    with open(cov_script, "w") as f:
        f.write("""\
import sys, glob, csv
from pathlib import Path
import numpy as np

ref = sys.argv[1]
out = sys.argv[2]

def parse_fasta_lengths(path):
    lengths, hid, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                hid = line[1:].split()[0]
                lengths[hid] = 0
            elif hid is not None:
                lengths[hid] += len(line.strip())
    return lengths

lengths = parse_fasta_lengths(ref)
# depth[insert][pool] = int32 per-bp array
depth = {c: {} for c in lengths}

for bg in sorted(glob.glob("cov__*.bedgraph")):
    pool = Path(bg).name[len("cov__"): -len(".bedgraph")]
    with open(bg) as fh:
        for line in fh:
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) != 4:
                continue
            c, s, e, v = parts
            if c not in lengths:
                continue
            s, e, v = int(s), int(e), int(v)
            arr = depth[c].get(pool)
            if arr is None:
                arr = np.zeros(lengths[c], dtype=np.int32)
                depth[c][pool] = arr
            arr[s:e] += v

def stats(arr):
    L = len(arr)
    return dict(
        length=L,
        mapped_bp=int(arr.sum()),
        mean_depth=round(float(arr.mean()), 4) if L else 0.0,
        median_depth=float(np.median(arr)) if L else 0.0,
        breadth_1x=round(float((arr >= 1).mean()), 6) if L else 0.0,
        breadth_5x=round(float((arr >= 5).mean()), 6) if L else 0.0,
    )

rows = []
for c in sorted(lengths):
    total = np.zeros(lengths[c], dtype=np.int64)
    for pool, arr in depth[c].items():
        total += arr
        st = stats(arr)
        rows.append([c, pool, st["length"], st["mapped_bp"], st["mean_depth"],
                     st["median_depth"], st["breadth_1x"], st["breadth_5x"]])
    st = stats(total.astype(np.int32) if total.max() < 2**31 else total)
    # recompute breadth on int64 to avoid overflow edge-cases
    L = lengths[c]
    rows.append([c, "all", L, int(total.sum()),
                 round(float(total.mean()), 4) if L else 0.0,
                 float(np.median(total)) if L else 0.0,
                 round(float((total >= 1).mean()), 6) if L else 0.0,
                 round(float((total >= 5).mean()), 6) if L else 0.0])

with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["insert_id", "pool", "length", "mapped_bp", "mean_depth",
                "median_depth", "breadth_1x", "breadth_5x"])
    w.writerows(rows)
n_covered = sum(1 for c in depth if depth[c])
print(f"inserts={len(lengths)} with_coverage={n_covered} rows={len(rows)}")
""")
    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {cov_script} {iref.container} {ocov.container}",
    )

    return ExecutionResult(
        manifest=[{out_bam: obam.local, out_csi: ocsi.local, out_cov: ocov.local}],
        success=obam.local.exists() and ocov.local.exists(),
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
