from metasmith.python_api import *
import json

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::hifiasm.oci"))
reads   = model.AddRequirement(lib.GetType("sequences::100x_long_reads"))
out     = model.AddProduct(lib.GetType("sequences::isolate_assembly"))

def protocol(context: ExecutionContext):
    ireads = context.Input(reads)
    iout = context.Output(out)

    threads = context.params.get('cpus')
    threads = "" if threads is None else f"-t{threads}"
    # Assemble inbred/homozygous genomes (-l0 disables duplication purging)
    context.ExecWithContainer(
        image = image,
        cmd = f"""
        hifiasm -o {iout.container} {threads} -l0 {ireads.container}
        """
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
    group_by=reads,
    resources=Resources(
        cpus=8,
        memory=Size.GB(32),
        duration=Duration(hours=18),
    )
)
