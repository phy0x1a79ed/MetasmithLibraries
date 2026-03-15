from metasmith.python_api import *

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("containers::eggnog-mapper.oci"))
proteins = model.AddRequirement(lib.GetType("transcriptomics::braker3_proteins"))
db = model.AddRequirement(lib.GetType("annotation::eggnog_db"))
out_results = model.AddProduct(lib.GetType("annotation::eggnog_results"))


def protocol(context: ExecutionContext):
    iproteins = context.Input(proteins)
    idb = context.Input(db)
    iout = context.Output(out_results)

    cpus = context.params.get("cpus")
    cpus = 16 if cpus is None else cpus

    # Strip terminal stop codons to keep emapper parsing robust.
    context.ExecWithContainer(
        image=image,
        binds=[(idb.external, "/db")],
        cmd=f"""
            sed '/^[^>]/s/\\*//g' {iproteins.container} > proteins.clean.faa
            emapper.py \
                --data_dir /db \
                -i proteins.clean.faa \
                -o results \
                --cpu {cpus} \
                -m diamond \
                --override
        """,
    )

    context.LocalShell(f"cp results.emapper.annotations {iout.local}")

    return ExecutionResult(
        manifest=[{out_results: iout.local}],
        success=iout.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=proteins,
    resources=Resources(
        cpus=16,
        memory=Size.GB(32),
        duration=Duration(hours=24),
    ),
)
