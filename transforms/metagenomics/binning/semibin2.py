from metasmith.python_api import *

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("containers::semibin.oci"))
asm         = model.AddRequirement(lib.GetType("sequences::assembly"))
bam         = model.AddRequirement(lib.GetType("alignment::bam"))
bins        = model.AddProduct(lib.GetType("binning::bin_directory"))
table       = model.AddProduct(lib.GetType("binning::contig_to_bin_table"))

def protocol(context: ExecutionContext):
    iasm = context.Input(asm)
    ibam = context.Input(bam)
    obins = context.Output(bins)
    otable = context.Output(table)

    threads = context.params.get('cpus', 8)
    workdir = "semibin_out"

    # Use global environment model (works for most samples)
    environment = "global"

    context.ExecWithContainer(
        image = image,
        cmd = f"""
            SemiBin2 single_easy_bin \
                -i {iasm.container} \
                -b {ibam.container} \
                -o {workdir} \
                --environment {environment} \
                -t {threads}
            mv {workdir}/output_bins {obins.container}
            cp {workdir}/contig_bins.tsv {otable.container}
        """
    )

    return ExecutionResult(
        manifest=[
            {
                bins: obins.local,
                table: otable.local,
            },
        ],
        success=obins.local.exists() and otable.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=8,
        memory=Size.GB(16),
        duration=Duration(hours=4),
    )
)
