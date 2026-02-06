from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::kofamscan.oci"))
ref     = model.AddRequirement(lib.GetType("functional_annotation::kofam_db"))
orfs    = model.AddRequirement(lib.GetType("sequences::orfs"))
ko      = model.AddProduct(lib.GetType("functional_annotation::kegg_orthologs"))

def protocol(context: ExecutionContext):
    iref  = context.Input(ref)
    iorfs = context.Input(orfs)
    iko   = context.Output(ko)

    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--cpu {threads}"

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            exec_annotation \
                -p {iref.container}/profiles \
                -k {iref.container}/ko_list \
                -o {iko.container} \
                --format detail-tsv \
                {threads} \
                {iorfs.container}
        """,
    )

    return ExecutionResult(
        manifest=[
            {
                ko: iko.local,
            },
        ],
        success=iko.local.exists(),
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
