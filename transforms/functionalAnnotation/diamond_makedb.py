from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::diamond.oci"))
fasta   = model.AddRequirement(lib.GetType("functional_annotation::protein_ref_fasta"))
db      = model.AddProduct(lib.GetType("functional_annotation::diamond_db"))

def protocol(context: ExecutionContext):
    ifasta = context.Input(fasta)
    idb    = context.Output(db)

    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--threads {threads}"

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            diamond makedb \
                --in {ifasta.container} \
                --db {idb.container} \
                {threads}
        """,
    )

    return ExecutionResult(
        manifest=[
            {
                db: idb.local,
            },
        ],
        success=idb.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=fasta,
    resources=Resources(
        cpus=4,
        memory=Size.GB(8),
        duration=Duration(hours=2),
    )
)
