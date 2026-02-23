from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
pair    = model.AddRequirement(lib.GetType("sequences::read_pair"))
r1      = model.AddRequirement(lib.GetType("sequences::forward_ora_reads"), parents={pair})
r2      = model.AddRequirement(lib.GetType("sequences::reverse_ora_reads"), parents={pair})
orad    = model.AddRequirement(lib.GetType("containers::orad.oci"))
bbtools = model.AddRequirement(lib.GetType("containers::bbtools.oci"))
out     = model.AddProduct(lib.GetType("sequences::short_reads"))

def protocol(context: ExecutionContext):
    ir1 = context.Input(r1)
    ir2 = context.Input(r2)
    iout = context.Output(out)
    threads = context.params.get('cpus')
    threads = "" if threads is None else f"-p {threads}"

    # Decompress ORA files to FASTQ using orad
    context.ExecWithContainer(
        image=orad,
        cmd=f'''
        orad "{ir1.container}" "{ir2.container}"
        '''
    )

    # orad writes .fastq files alongside the .fastq.ora files
    # Derive the decompressed FASTQ paths by stripping the .ora extension
    r1_fq = f'"{ir1.container}"'.replace('.ora"', '"')
    r2_fq = f'"{ir2.container}"'.replace('.ora"', '"')

    # Interleave R1+R2 with reformat.sh and compress with pigz
    context.ExecWithContainer(
        image=bbtools,
        cmd=f'''
        reformat.sh \
            in1={r1_fq} \
            in2={r2_fq} \
            out=stdout.fq \
        | pigz {threads} > {iout.container}
        '''
    )

    return ExecutionResult(
        manifest=[{
            out: iout.local
        }],
        success=iout.local.exists()
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=r1,
    resources=Resources(
        cpus=2,
        memory=Size.GB(8),
        duration=Duration(hours=3),
    )
)
