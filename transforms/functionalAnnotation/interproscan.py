from metasmith.python_api import *
from pathlib import Path

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("containers::interproscan.oci"))
orfs = model.AddRequirement(lib.GetType("sequences::orfs"))
data_dir = model.AddRequirement(lib.GetType("annotation::interproscan_data"))
out_json = model.AddProduct(lib.GetType("annotation::interproscan_json"))
out_gff = model.AddProduct(lib.GetType("annotation::interproscan_gff"))


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    idata = context.Input(data_dir)
    ijson = context.Output(out_json)
    igff = context.Output(out_gff)

    cpus = context.params.get("cpus")
    cpus_string = "" if cpus is None else f"--cpu {cpus}"

    # sed to remove stop codon asterisks from protein sequences
    # Run InterProScan with external data directory
    context.ExecWithContainer(
        image=image,
        binds=[(idata.external/"data", "/opt/interproscan/data")],
        cmd=f"""
            mkdir -p output
            sed '/^[^>]/s/\*//g' {iorfs.container} > ./{iorfs.container.stem}.clean.faa
            /opt/interproscan/interproscan.sh \
                --disable-precalc \
                --verbose \
                --seqtype p \
                {cpus_string} \
                -i ./{iorfs.container.stem}.clean.faa \
                -f json,gff3 \
                -d output
        """,
    )

    # Copy outputs
    json_files = list(Path("output").glob("*.json"))
    gff_files = list(Path("output").glob("*.gff3"))

    if json_files:
        context.LocalShell(f"cp {json_files[0]} {ijson.local}")
    if gff_files:
        context.LocalShell(f"cp {gff_files[0]} {igff.local}")

    return ExecutionResult(
        manifest=[
            {
                out_json: ijson.local,
                out_gff: igff.local,
            },
        ],
        success=ijson.local.exists() and igff.local.exists(),
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
