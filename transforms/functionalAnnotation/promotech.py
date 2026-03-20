from metasmith.python_api import *

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image    = model.AddRequirement(lib.GetType("containers::promotech.oci"))
assembly = model.AddRequirement(lib.GetType("sequences::assembly"))

out_pred = model.AddProduct(lib.GetType("functional_annotation::promotech_predictions"))


def protocol(context: ExecutionContext):
    iasm = context.Input(assembly)
    opred = context.Output(out_pred)

    # Step 1: Parse genome into 40nt sliding windows and encode
    # PromoTech uses relative model paths, so CWD must be /opt/promotech.
    # iasm.container is relative to /ws, so use absolute path.
    # Must cd back to /ws so the bounce script can write exitcode.
    context.ExecWithContainer(
        image=image,
        cmd=f"""
            cd /opt/promotech &&
            python promotech.py -pg -f /ws/{iasm.container} -o /tmp/pt -m RF-HOT &&
            cd /ws
        """,
    )

    # Step 2: Predict promoters on forward + inverse strands
    context.ExecWithContainer(
        image=image,
        cmd=f"""
            cd /opt/promotech &&
            python promotech.py -g -i /tmp/pt -o /tmp/pt -m RF-HOT -t 0.5 &&
            cd /ws
        """,
    )

    # Copy result TSV
    context.LocalShell(f"cp /tmp/pt/genome_predictions.csv {opred.local}")

    return ExecutionResult(
        manifest=[{out_pred: opred.local}],
        success=opred.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=assembly,
    resources=Resources(
        cpus=2,
        memory=Size.GB(48),
        duration=Duration(hours=4),
    ),
)
