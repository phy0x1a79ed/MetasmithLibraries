from metasmith.python_api import *
import json

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
img_mm2   = model.AddRequirement(lib.GetType("containers::minimap2.oci"))
img_sam   = model.AddRequirement(lib.GetType("containers::samtools.oci"))
img_vs    = model.AddRequirement(lib.GetType("containers::vsearch.oci"))
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
meta      = model.AddRequirement(lib.GetType("sequences::read_metadata"))
reads     = model.AddRequirement(lib.GetType("sequences::clean_short_reads"), parents={meta})
backbone  = model.AddRequirement(lib.GetType("fosmids::vector_backbone"))
out       = model.AddProduct(lib.GetType("fosmids::pool_size_estimate"))

def protocol(context: ExecutionContext):
    ireads = context.Input(reads)
    imeta  = context.Input(meta)
    iback  = context.Input(backbone)
    iout   = context.Output(out)

    with open(imeta.local) as j:
        read_meta = json.load(j)
    parity = read_meta["parity"]
    assert parity == "paired", f"estimate_pool_size requires paired reads, got [{parity}]"

    threads     = context.params.get("cpus")
    threads_mm2 = "" if threads is None else f"-t {threads}"
    threads_sam = "" if threads is None else f"-@ {threads}"
    threads_vs  = "" if threads is None else f"--threads {threads}"

    # Stage 1: prepare normalized backbone fasta + record its 7-nt signature
    norm_script = "norm_backbone.py"
    with open(norm_script, "w") as f:
        f.write("""\
import sys, json
from Bio import SeqIO
src = sys.argv[1]
dst = sys.argv[2]
sig_path = sys.argv[3]

rec = next(SeqIO.parse(src, "fasta"))
seq = str(rec.seq)
sig = seq[:7]
with open(dst, "w") as out:
    out.write(">vector_backbone\\n")
    s = seq
    for i in range(0, len(s), 80):
        out.write(s[i:i+80] + "\\n")
with open(sig_path, "w") as j:
    json.dump({"signature": sig}, j)
print(f"signature={sig}")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {norm_script} {iback.container} backbone.fasta signature.json",
    )

    # Stage 2: map reads to backbone, keep mapped only as fastq
    context.ExecWithContainer(
        image=img_mm2,
        cmd=f"""\
            minimap2 -a -x sr {threads_mm2} --secondary=no \
                backbone.fasta {ireads.container} > mapped.sam
        """,
    )

    context.ExecWithContainer(
        image=img_sam,
        cmd=f"""\
            samtools view -ub -F 4 {threads_sam} mapped.sam \
            | samtools fastq {threads_sam} -N - >mapped.fastq
            rm -f mapped.sam
        """,
    )

    # Stage 3: trim signature, keep last 100bp, write candidates fasta
    trim_script = "trim_signatures.py"
    with open(trim_script, "w") as f:
        f.write("""\
import json
from Bio import SeqIO

with open("signature.json") as j:
    sig = json.load(j)["signature"]

CUT = 100
records = list(SeqIO.parse("mapped.fastq", "fastq"))

# keep both orientations
def _candidates():
    for i, rec in enumerate(records):
        seq = str(rec.seq)
        yield i, False, seq
        yield i, True, str(rec.seq.reverse_complement())

kept = []
with open("candidates.fasta", "w") as out:
    n = 0
    for i, is_rc, seq in _candidates():
        if sig not in seq: continue
        new_seq = seq[: seq.index(sig)]
        if len(new_seq) < CUT: continue
        trim = new_seq[len(new_seq) - CUT: len(new_seq)]
        out.write(f">{n}\\n{trim}\\n")
        kept.append((i, is_rc))
        n += 1

with open("kept_index.json", "w") as j:
    json.dump(kept, j)
print(f"candidates={len(kept)}")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {trim_script}",
    )

    # Stage 4: vsearch cluster at 90% identity
    context.ExecWithContainer(
        image=img_vs,
        cmd=f"""\
            if [ -s candidates.fasta ]; then
                vsearch -sizeout -id 0.9 {threads_vs} \
                    -cluster_fast candidates.fasta \
                    -uc clusters.uc \
                    -consout clusters.fasta \
                    -centroids centroids.fasta >vsearch.log 2>&1
            else
                : >centroids.fasta
            fi
        """,
    )

    # Stage 5: count clusters, write estimate JSON
    count_script = "count_clusters.py"
    with open(count_script, "w") as f:
        f.write("""\
import json, sys
from Bio import SeqIO

out_path = sys.argv[1]

estimate = 0
size_ones = 0
try:
    for rec in SeqIO.parse("centroids.fasta", "fasta"):
        estimate += 1
        # vsearch encodes cluster size as ;size=<n>
        parts = rec.id.split(";")
        size_field = next((p for p in parts if p.startswith("size=")), "size=0")
        try:
            sz = int(size_field.split("=", 1)[1])
        except ValueError:
            sz = 0
        if sz == 1:
            size_ones += 1
except FileNotFoundError:
    pass

total = estimate + size_ones
result = dict(size=estimate, size_with_singletons=total)
with open(out_path, "w") as j:
    json.dump(result, j)
print(json.dumps(result))
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {count_script} {iout.container}",
    )

    return ExecutionResult(
        manifest=[
            {
                out: iout.local,
            },
        ],
        success=iout.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=meta,
    resources=Resources(
        cpus=4,
        memory=Size.GB(16),
        duration=Duration(hours=4),
    )
)
