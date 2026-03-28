from metasmith.python_api import *
from pathlib import Path
import shutil

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("containers::interproscan.oci"))
orfs = model.AddRequirement(lib.GetType("sequences::orfs"))
data_dir = model.AddRequirement(lib.GetType("annotation::interproscan_data"))
out_gff = model.AddProduct(lib.GetType("annotation::interproscan_gff3"))
out_json = model.AddProduct(lib.GetType("annotation::interproscan_json"))


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    idata = context.Input(data_dir)
    igff = context.Output(out_gff)
    ijson = context.Output(out_json)

    cpus = context.params.get("cpus")
    cpus_string = "" if cpus is None else f"--cpu {cpus}"

    context.ExecWithContainer(
        image=image,
        binds=[(idata.external/"data", "/opt/interproscan/data")],
        cmd=f"""
            mkdir -p output
            sed '/^[^>]/s/\\*//g' {iorfs.container} > ./{iorfs.container.stem}.clean.faa
            /opt/interproscan/interproscan.sh \
                --disable-precalc \
                --seqtype p \
                --appl Pfam,NCBIfam,HAMAP,SFLD,MobiDBLite \
                {cpus_string} \
                -i ./{iorfs.container.stem}.clean.faa \
                -f gff3,json \
                -d output \
                -dp
        """,
    )

    gff_files = list(Path("output").glob("*.gff3"))
    if gff_files:
        shutil.copy2(str(gff_files[0]), str(igff.local))

    json_files = list(Path("output").glob("*.json"))
    if json_files:
        shutil.copy2(str(json_files[0]), str(ijson.local))

    return ExecutionResult(
        manifest=[
            {
                out_gff: igff.local,
                out_json: ijson.local,
            },
        ],
        success=igff.local.exists() and ijson.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=8,
        memory=Size.GB(16),
        duration=Duration(hours=24),
    ),
)
