"""Per-pool read depth over the clustered contigs -> a merged coverage-profile artifact.

Maps each pool's reads (staged one file per pool, named `<pool>.fq`) to the
clustered contigs, converts each alignment to a per-bp bedgraph, then folds all
pools into one `fosmids::pool_coverage_profiles` pickle:
    {contig -> {pool -> int32 depth array of length L}}
keyed on the original clustered contig ids.  This mirrors scadc
`coverage_profile.build_profiles` (a pool's barcodes are summed upstream, so the
per-pool track here is that pool's summed depth).  The profile is immutable
downstream: chimera_split / coverage_trim slice it by each contig's header interval.

Aggregation: one `fosmids::recovery_experiment` node groups the run; the clustered
contigs and every pool's reads are parented to it, so `group_by=exp` gives ONE job
the single clustered fasta plus all pools' reads via `InputGroup`.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
clustered = model.AddRequirement(lib.GetType("fosmids::clustered_contigs"), parents={exp})
reads     = model.AddRequirement(lib.GetType("sequences::clean_short_reads"), parents={exp})
img_mm2   = model.AddRequirement(lib.GetType("containers::minimap2.oci"))
img_sam   = model.AddRequirement(lib.GetType("containers::samtools.oci"))
img_bed   = model.AddRequirement(lib.GetType("containers::bedtools.oci"))
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
out_prof  = model.AddProduct(lib.GetType("fosmids::pool_coverage_profiles"))

def protocol(context: ExecutionContext):
    iclustered  = context.Input(clustered)
    reads_paths = context.InputGroup(reads)
    oprof       = context.Output(out_prof)

    threads     = context.params.get("cpus")
    threads_mm2 = "" if threads is None else f"-t {threads}"
    threads_sam = "" if threads is None else f"-@ {threads}"

    # one mapping -> bedgraph per staged pool file (pool id = filename stem)
    for rp in reads_paths:
        pool = rp.local.name
        for suf in ("_unmapped", ".fq", ".fastq", ".gz", ".fasta", ".fna"):
            pool = pool.replace(suf, "")
        context.ExecWithContainer(
            image=img_mm2,
            cmd=f"""\
                minimap2 -a -x sr {threads_mm2} --secondary=no \
                    {iclustered.container} {rp.container} > aln.{pool}.sam
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
                rm -f aln.{pool}.bam
            """,
        )

    # fold every cov__<pool>.bedgraph into per-contig x per-pool per-bp arrays
    fold_script = "fold_profiles.py"
    with open(fold_script, "w") as f:
        f.write(f"""\
import sys, glob, gzip, pickle
from pathlib import Path
import numpy as np
from Bio import SeqIO

clustered = sys.argv[1]
out_path  = sys.argv[2]

lengths = {{s.id: len(s) for s in SeqIO.parse(clustered, "fasta")}}
profiles = {{c: {{}} for c in lengths}}

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
            arr = profiles[c].get(pool)
            if arr is None:
                arr = np.zeros(lengths[c], dtype=np.int32)
                profiles[c][pool] = arr
            arr[s:e] += v   # '+=' so repeated pool files accumulate

with gzip.open(out_path, "wb") as gz:
    pickle.dump(dict(profiles=profiles, lengths=lengths), gz,
                protocol=pickle.HIGHEST_PROTOCOL)
n_cov = sum(1 for c in profiles if profiles[c])
print(f"contigs={{len(lengths)}} with_coverage={{n_cov}}")
""")
    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {fold_script} {iclustered.container} {oprof.container}",
    )

    return ExecutionResult(
        manifest=[{out_prof: oprof.local}],
        success=oprof.local.exists(),
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
