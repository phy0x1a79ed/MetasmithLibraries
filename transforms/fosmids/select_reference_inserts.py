"""Select the canonical fosmid reference: putative inserts with length >= cutoff.

The recover-fosmids DAG emits 199 `putative_inserts`; every downstream analysis
(read-mapping, 4-lane annotation, ECSPr) runs on the subset with length >=
`min_insert_length` (default 29000 bp -- the phage-packaging / lambda-window cutoff
that matches scadc's `fabfos_map` input of 132 contigs).  This transform makes that
selection reproducible in-DAG rather than as a one-off hand filter.

Filter on sequence LENGTH, not the report's `band` string (the report's trailing
column can carry CRLF, which breaks naive `band != short` matching).
"""
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
inserts  = model.AddRequirement(lib.GetType("fosmids::putative_inserts"), parents={exp})
img_pyds = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
out_ref  = model.AddProduct(lib.GetType("fosmids::reference_inserts"))

def protocol(context: ExecutionContext):
    iins = context.Input(inserts)
    oref = context.Output(out_ref)

    min_len = context.params.get("min_insert_length", 29000)

    select_script = "select_reference_inserts.py"
    with open(select_script, "w") as f:
        f.write(f"""\
import sys

MIN_LEN = {min_len}
src, out = sys.argv[1], sys.argv[2]

def parse_fasta(path):
    hdr, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line, []
            else:
                seq.append(line)
    if hdr is not None:
        yield hdr, "".join(seq)

kept = 0
with open(out, "w") as w:
    for hdr, seq in parse_fasta(src):
        if len(seq) >= MIN_LEN:
            w.write(hdr + "\\n")
            for i in range(0, len(seq), 60):
                w.write(seq[i:i+60] + "\\n")
            kept += 1
print(f"reference_inserts: {{kept}} records with length >= {{MIN_LEN}}")
""")
    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {select_script} {iins.container} {oref.container}",
    )

    return ExecutionResult(
        manifest=[{out_ref: oref.local}],
        success=oref.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=inserts,
    resources=Resources(
        cpus=1,
        memory=Size.GB(2),
        duration=Duration(minutes=10),
    )
)
