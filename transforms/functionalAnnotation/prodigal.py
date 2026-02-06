from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::prodigal.oci"))
asm     = model.AddRequirement(lib.GetType("sequences::assembly"))
orfs    = model.AddProduct(lib.GetType("sequences::orfs"))
gff     = model.AddProduct(lib.GetType("sequences::gff"))

def protocol(context: ExecutionContext):
    iasm  = context.Input(asm)
    iorfs = context.Output(orfs)
    igff  = context.Output(gff)

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            prodigal \
                -i {iasm.container} \
                -a {iorfs.container} \
                -o {igff.container} \
                -f gff \
                -p meta
        """,
    )

    return ExecutionResult(
        manifest=[
            {
                orfs: iorfs.local,
                gff: igff.local,
            },
        ],
        success=iorfs.local.exists() and igff.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=1,
        memory=Size.GB(2),
        duration=Duration(hours=1),
    )
)
