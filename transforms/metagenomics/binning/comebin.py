from metasmith.python_api import *

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("containers::comebin.oci"))
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
    workdir = "comebin_out"
    bam_dir = "bam_input"

    context.ExecWithContainer(
        image = image,
        cmd = f"""
            mkdir -p {bam_dir}
            cp -L {ibam.container} {bam_dir}/
            run_comebin.sh -a {iasm.container} -o {workdir} -p {bam_dir} -t {threads}
            mv {workdir}/comebin_res/comebin_res_bins {obins.container}
            cp {workdir}/comebin_res/comebin_res_contig_bin.tsv {otable.container}
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
        memory=Size.GB(32),
        duration=Duration(hours=12),
    )
)
