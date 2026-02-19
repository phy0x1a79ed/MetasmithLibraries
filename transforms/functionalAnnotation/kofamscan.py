from metasmith.python_api import *

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image    = model.AddRequirement(lib.GetType("containers::kofamscan.oci"))
orfs     = model.AddRequirement(lib.GetType("sequences::orfs"))
profiles = model.AddRequirement(lib.GetType("annotation::kofamscan_profiles"))
ko_list  = model.AddRequirement(lib.GetType("annotation::kofamscan_ko_list"))
out_results = model.AddProduct(lib.GetType("annotation::kofamscan_results"))


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    iprofiles = context.Input(profiles)
    iko_list = context.Input(ko_list)
    iout = context.Output(out_results)

    threads = context.params.get("cpus", 8)

    # Run KofamScan with external database references
    context.ExecWithContainer(
        image=image,
        binds=[
            (iprofiles.external, "/profiles"),
            (iko_list.external.parent, "/ko"),
        ],
        cmd=f"""
            exec_annotation \
                -o kofam_results.txt \
                --profile=/profiles \
                --ko-list=/ko/{iko_list.external.name} \
                --cpu={threads} \
                --e-value=0.01 \
                --format=detail \
                --no-report-unannotated \
                {iorfs.container}
        """,
    )

    # Copy output
    context.LocalShell(f"cp kofam_results.txt {iout.local}")

    return ExecutionResult(
        manifest=[
            {
                out_results: iout.local,
            },
        ],
        success=iout.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=8,
        memory=Size.GB(16),
        duration=Duration(hours=8),
    ),
)
