from metasmith.python_api import *
import json

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
img_mm2 = model.AddRequirement(lib.GetType("containers::minimap2.oci"))
img_sam = model.AddRequirement(lib.GetType("containers::samtools.oci"))
meta    = model.AddRequirement(lib.GetType("sequences::read_metadata"))
reads   = model.AddRequirement(lib.GetType("sequences::clean_short_reads"), parents={meta})
host    = model.AddRequirement(lib.GetType("sequences::background_genome"))
out     = model.AddProduct(lib.GetType("sequences::host_filtered_short_reads"))

def protocol(context: ExecutionContext):
    ireads = context.Input(reads)
    imeta  = context.Input(meta)
    ihost  = context.Input(host)
    iout   = context.Output(out)

    with open(imeta.local) as j:
        read_meta = json.load(j)
    parity = read_meta["parity"]
    assert parity in {"single", "paired"}, f"unknown parity: [{parity}]"
    sr_params = "-x sr" if parity == "paired" else ""

    threads = context.params.get("cpus")
    threads_mm2 = "" if threads is None else f"-t {threads}"
    threads_sam = "" if threads is None else f"-@ {threads}"

    bam = "temp.bam"
    context.ExecWithContainer(
        image=img_mm2,
        cmd=f"""\
            minimap2 -a {sr_params} {threads_mm2} --secondary=no \
                {ihost.container} {ireads.container} > temp.sam
        """,
    )
    context.ExecWithContainer(
        image=img_sam,
        cmd=f"""\
            samtools view -ub -f 4 {threads_sam} temp.sam \
            | samtools fastq {threads_sam} -N - \
            | gzip > {iout.container}
            rm -f temp.sam {bam}
        """,
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
    group_by=meta,
    resources=Resources(
        cpus=4,
        memory=Size.GB(16),
        duration=Duration(hours=12),
    )
)
