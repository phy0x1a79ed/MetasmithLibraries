from metasmith.python_api import *
import json

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::megahit.oci"))
meta    = model.AddRequirement(lib.GetType("sequences::read_metadata"))
reads   = model.AddRequirement(lib.GetType("sequences::clean_short_reads"), parents={meta})
out     = model.AddProduct(lib.GetType("sequences::assembly"))

def protocol(context: ExecutionContext):
    ireads=context.Input(reads)
    imeta=context.Input(meta)
    iout=context.Output(out)
    with open(imeta.local) as j:
        read_meta = json.load(j)
    parity = read_meta["parity"]
    assert parity in {"single", "paired"}, f"unknown parity: [{parity}]"
    if parity == "paired":
        parg = "--12"
    else:
        parg = "-r"
    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--num-cpu-threads {threads}"
    # memory defaults to 0.9 of available [--memory 0.9]
    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            megahit {threads} \
                {parg} {ireads.container} \
                -o megahit_ws
            [[ $(head megahit_ws/final.contigs.fa | wc --chars) -ne 0 ]] && mv megahit_ws/final.contigs.fa {iout.container} || echo "assembly was empty"
        """,
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
        cpus=8,
        memory=Size.GB(32),
        duration=Duration(hours=18),
    )
)
