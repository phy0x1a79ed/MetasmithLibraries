from metasmith.python_api import *

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("containers::eggnog-mapper.oci"))
orfs = model.AddRequirement(lib.GetType("sequences::orfs"))
db = model.AddRequirement(lib.GetType("annotation::eggnog_db"))
out_results = model.AddProduct(lib.GetType("annotation::eggnog_mapper_results"))

def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    idb = context.Input(db)
    iout = context.Output(out_results)

    threads = context.params.get("cpus", 8)

    context.ExecWithContainer(
        image=image,
        binds=[(idb.external, "/db")],
        cmd=f"""
            emapper.py \
                -i {iorfs.container} \
                --output eggnog_out \
                --cpu {threads} \
                --data_dir /db \
                -m diamond \
                --override
        """,
    )

    # eggnog-mapper outputs eggnog_out.emapper.annotations
    context.LocalShell(f"cp eggnog_out.emapper.annotations {iout.local}")

    return ExecutionResult(
        manifest=[{out_results: iout.local}],
        success=iout.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=8,
        memory=Size.GB(32),
        duration=Duration(hours=12),
    ),
)
