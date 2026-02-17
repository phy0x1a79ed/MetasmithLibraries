from metasmith.python_api import *

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("containers::metabat2.oci"))
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
    depth_file = "depth.txt"
    bin_prefix = "metabat_bins/bin"

    context.ExecWithContainer(
        image = image,
        cmd = f"""
            mkdir -p metabat_bins
            jgi_summarize_bam_contig_depths --outputDepth {depth_file} {ibam.container}
            metabat2 -i {iasm.container} -a {depth_file} -o {bin_prefix} -t {threads}
            mv metabat_bins {obins.container}

            # Generate contig-to-bin table from bin files
            echo -e "contig\\tbin" > {otable.container}
            for f in {obins.container}/*.fa; do
                bin_name=$(basename "$f" .fa)
                grep "^>" "$f" | sed 's/^>//' | while read contig; do
                    echo -e "$contig\\t$bin_name"
                done
            done >> {otable.container}
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
        memory=Size.GB(8),
        duration=Duration(hours=2),
    )
)
