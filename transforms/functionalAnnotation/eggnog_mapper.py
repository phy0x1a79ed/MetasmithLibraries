from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
image    = model.AddRequirement(lib.GetType("containers::eggnog-mapper.oci"))
orfs     = model.AddRequirement(lib.GetType("sequences::orfs"))
data_dir = model.AddRequirement(lib.GetType("annotation::eggnog_data"))
out      = model.AddProduct(lib.GetType("annotation::eggnog_results"))

def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    idata = context.Input(data_dir)
    iout  = context.Output(out)

    cpus = context.params.get("cpus")
    cpus = 16 if cpus is None else cpus

    # Strip stop codon asterisks from protein sequences (same as interproscan)
    context.ExecWithContainer(
        image=image,
        binds=[(idata.external, "/eggnog_data")],
        cmd=f"""
            sed '/^[^>]/s/\\*//g' {iorfs.container} > ./proteins.clean.faa
            emapper.py \
                --data_dir /eggnog_data \
                -i ./proteins.clean.faa \
                -o results \
                -m diamond \
                --cpu {cpus}
        """,
    )

    context.LocalShell(f"cp results.emapper.annotations {iout.local}")

    return ExecutionResult(
        manifest=[{out: iout.local}],
        success=iout.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=16,
        memory=Size.GB(32),
        duration=Duration(hours=24),
    ),
)
