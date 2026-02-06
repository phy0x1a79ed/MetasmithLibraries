from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::diamond.oci"))
ref     = model.AddRequirement(lib.GetType("functional_annotation::diamond_db"))
orfs    = model.AddRequirement(lib.GetType("sequences::orfs"))
hits    = model.AddProduct(lib.GetType("functional_annotation::diamond_hits"))

def protocol(context: ExecutionContext):
    iref  = context.Input(ref)
    iorfs = context.Input(orfs)
    ihits = context.Output(hits)

    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--threads {threads}"

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            diamond blastp \
                --query {iorfs.container} \
                --db {iref.container} \
                --out {ihits.container} \
                --outfmt 6 \
                --sensitive \
                {threads}
        """,
    )

    return ExecutionResult(
        manifest=[
            {
                hits: ihits.local,
            },
        ],
        success=ihits.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=4,
        memory=Size.GB(8),
        duration=Duration(hours=2),
    )
)
