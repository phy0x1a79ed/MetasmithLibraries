from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
img_blast = model.AddRequirement(lib.GetType("containers::blast.oci"))
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
asm       = model.AddRequirement(lib.GetType("sequences::assembly"))
kept      = model.AddProduct(lib.GetType("sequences::length_filtered_assembly"))
discarded = model.AddProduct(lib.GetType("sequences::discarded_contigs"))

def protocol(context: ExecutionContext):
    iasm  = context.Input(asm)
    ikept = context.Output(kept)
    idisc = context.Output(discarded)

    min_length = context.params.get("min_length", 1000)
    max_pident = context.params.get("max_pident", 99)
    threads    = context.params.get("cpus")

    threads_arg = "" if threads is None else f"-num_threads {threads}"

    # Step 1: length filter into a temp fasta and build blast db
    prefilter_script = "prefilter.py"
    with open(prefilter_script, "w") as f:
        f.write(f"""\
import sys
from Bio import SeqIO

src = sys.argv[1]
dst = sys.argv[2]
min_len = {min_length}

with open(dst, "w") as out:
    for rec in SeqIO.parse(src, "fasta"):
        if len(rec.seq) < min_len:
            continue
        out.write(f">{{rec.id}} length={{len(rec.seq)}}\\n")
        s = str(rec.seq)
        for i in range(0, len(s), 80):
            out.write(s[i:i+80] + "\\n")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {prefilter_script} {iasm.container} prefiltered.fna",
    )

    context.ExecWithContainer(
        image=img_blast,
        cmd=f"""\
            mkdir -p blast_db
            makeblastdb -dbtype nucl -in prefiltered.fna -out blast_db/contigs >makeblastdb.log 2>&1
            blastn -perc_identity 50 {threads_arg} \
                -query prefiltered.fna \
                -db blast_db/contigs \
                -outfmt "6 qseqid sseqid nident length qlen slen qstart qend sstart send" \
                -out pairwise.tsv >blastn.log 2>&1
        """,
    )

    # Step 2: parse hits and split kept/discarded
    split_script = "split_dedupe.py"
    with open(split_script, "w") as f:
        f.write(f"""\
import sys
import pandas as pd
from Bio import SeqIO

prefiltered = sys.argv[1]
pairwise = sys.argv[2]
kept_path = sys.argv[3]
discarded_path = sys.argv[4]
max_pident = {max_pident}

cols = "qseqid sseqid nident length qlen slen qstart qend sstart send".split()
to_remove = set()
try:
    df = pd.read_csv(pairwise, sep="\\t", header=None, names=cols)
    similars = {{}}
    for _, row in df.iterrows():
        if row.qseqid == row.sseqid:
            continue
        if row.qlen == 0:
            continue
        if (row.nident / row.qlen) * 100 >= max_pident:
            similars.setdefault(row.sseqid, set()).add(row.qseqid)
    for k, members in sorted(similars.items(), key=lambda t: len(t[1]), reverse=True):
        if k in to_remove:
            continue
        to_remove |= members
except (pd.errors.EmptyDataError, FileNotFoundError):
    pass

records = list(SeqIO.parse(prefiltered, "fasta"))
records.sort(key=lambda r: len(r.seq), reverse=True)
import math
places = max(1, int(math.log10(len(records) or 1)) + 1)

ki = di = 0
with open(kept_path, "w") as kf, open(discarded_path, "w") as df_out:
    for rec in records:
        if rec.id in to_remove:
            di += 1
            i = di
            target = df_out
        else:
            ki += 1
            i = ki
            target = kf
        target.write(f">C{{i:0{{places}}}} length={{len(rec.seq)}}\\n")
        s = str(rec.seq)
        for j in range(0, len(s), 80):
            target.write(s[j:j+80] + "\\n")

print(f"kept={{ki}} discarded={{di}}")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {split_script} prefiltered.fna pairwise.tsv {ikept.container} {idisc.container}",
    )

    return ExecutionResult(
        manifest=[
            {
                kept: ikept.local,
                discarded: idisc.local,
            },
        ],
        success=ikept.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=4,
        memory=Size.GB(8),
        duration=Duration(hours=4),
    )
)
