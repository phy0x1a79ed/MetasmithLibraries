from pathlib import Path
from metasmith.python_api import *
import re

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
pan     = model.AddRequirement(lib.GetType("pangenome::pangenome"))
gbk     = model.AddRequirement(lib.GetType("sequences::gbk"), parents={pan})
image   = model.AddRequirement(lib.GetType("env::ppanggolin.env"))
matrix  = model.AddProduct(lib.GetType("pangenome::ppanggolin_matrix"))
pg      = model.AddProduct(lib.GetType("pangenome::ppanggolin_raw"))

def protocol(context: ExecutionContext):
    dep_paths=context.InputGroup(gbk)

    gb_list = Path("genbank_manifest.list")
    with open(gb_list, "w") as f:
        for p in dep_paths:
            name = None
            with open(p.local) as g:
                for i, l in enumerate(g):
                    if i > 15: continue
                    if "DEFINITION" not in l: continue
                    for g in re.finditer(r"substr\.?\s?([^\s]+)|strain\s?([^\s]+)", l):
                        a, b = g.group(1), g.group(2)
                        name = a if a else b
                        if name: break
                    if not name: name = l.replace("DEFINITION", "")
                    name = name.strip()
                    name = name.replace(" ", "-") # the regex should take care of this already...
                    for x in ",.'\"":
                        name = name.replace(x, "")
                    break
            if not name: name = p.local.name
            f.write(name+"\t"+str(p.container)+"\n")

    ipg = context.Output(pg)
    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--cpu {threads}"
    context.ExecWithContainer(
        image=image,
        cmd=f"ppanggolin all --anno {gb_list} {threads} --output {ipg.container}",
    )
    imatrix = context.Output(matrix)
    context.LocalShell(f"cp {ipg.local}/matrix.csv {imatrix.local}")
    return ExecutionResult(
        manifest=[
            {
                pg: ipg.local,
                matrix: imatrix.local
            },
        ],
        success=ipg.local.exists() and imatrix.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=pan,
    resources=Resources(
        cpus=4,
        memory=Size.GB(8),
        duration=Duration(hours=3),
    )
)
