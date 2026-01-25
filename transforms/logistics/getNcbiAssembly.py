import glob
import os
from pathlib import Path
import shutil
from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
dep     = model.AddRequirement(lib.GetType("ncbi::assembly_accession"))
image   = model.AddRequirement(lib.GetType("containers::ncbi-datasets.oci"))
fna     = model.AddProduct(lib.GetType("sequences::assembly"))
faa     = model.AddProduct(lib.GetType("sequences::orfs"))
gff     = model.AddProduct(lib.GetType("sequences::gff"))
gbk     = model.AddProduct(lib.GetType("sequences::gbk"))

def protocol(context: ExecutionContext):
    dep_path=context.Input(dep)

    with open(dep_path.local) as f:
        acc = f.readline().strip()

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            datasets download genome accession {acc} \
                --include gff3,protein,genome,gbff
        """,
    )
    context.LocalShell(f"unzip ncbi_dataset.zip")

    output_manifest = {}
    def fix_out(dep, p: Path):
        op = context.Output(dep)
        shutil.move(p, op.local)
        output_manifest[dep] = op.local
    for f in glob.glob("ncbi_dataset/*/*/*"):
        p = Path(f)
        Log.Info(f"scanning file [{p}]")
        match(p.name):
            case "genomic.gff":
                fix_out(gff, p)
            case "genomic.gbff":
                fix_out(gbk, p)
            case "protein.faa":
                fix_out(faa, p)
        if not p.name.startswith("cds") and p.name.endswith("genomic.fna"):
                fix_out(fna, p)
    return ExecutionResult(
        manifest=[
            output_manifest,
        ],
        success=len(output_manifest)==len(model.produces[0]), # no branching
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=dep,
    labels=["local"],
    resources=Resources(
        cpus=1,
        memory=Size.GB(1),
    )
)
