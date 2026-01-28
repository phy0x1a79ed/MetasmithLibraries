from metasmith.python_api import *

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("containers::checkm.oci")) # database included at /checkm_database
asm         = model.AddRequirement(lib.GetType("sequences::assembly"))
out         = model.AddProduct(lib.GetType("taxonomy::checkm_stats"))
raw         = model.AddProduct(lib.GetType("taxonomy::checkm_raw"))

def protocol(context: ExecutionContext):
    iasm = context.Input(asm)
    iout = context.Output(out)
    iraw = context.Output(raw)

    threads = context.params.get('cpus')
    threads = "" if threads is None else f"-t {threads}"
    ext = iasm.container.suffix.replace(".", "")
    temp_ws = "checkm_ws"
    context.ExecWithContainer(
        image = image,
        cmd = f"""
            mkdir input
            cp -L {iasm.container} ./input/
            checkm lineage_wf {threads} -x {ext} ./input ./{temp_ws}
            checkm qa ./{temp_ws}/lineage.ms ./{temp_ws} --out_format 2 --tab_table --file {iout.container}
            mv ./{temp_ws} {iraw.container}
        """
    )
    
    return ExecutionResult(
        manifest=[
            {
                out: iout.local,
                raw: iraw.local,
            },
        ],
        success=iout.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=2,
        memory=Size.GB(44), # checkm needs 40
        duration=Duration(hours=1),
    )
)
