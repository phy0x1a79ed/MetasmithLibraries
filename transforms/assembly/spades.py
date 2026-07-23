from metasmith.python_api import *
import json

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("env::spades.env"))
meta    = model.AddRequirement(lib.GetType("sequences::read_metadata"))
reads   = model.AddRequirement(lib.GetType("sequences::clean_short_reads"), parents={meta})
out     = model.AddProduct(lib.GetType("sequences::spades_assembly"))

def protocol(context: ExecutionContext):
    ireads=context.Input(reads)
    imeta=context.Input(meta)
    iout=context.Output(out)
    with open(imeta.local) as j:
        read_meta = json.load(j)
    parity = read_meta["parity"]
    assert parity in {"single", "paired"}, f"unknown parity: [{parity}]"
    if parity == "paired":
        # metaSPAdes (--meta) is pooled-clone/metagenome-appropriate, but only
        # accepts paired input; single-end falls back to plain spades.py.
        mode = "--meta"
        reads_arg = f"--12 {ireads.container}"
    else:
        mode = ""
        reads_arg = f"-s {ireads.container}"

    threads = context.params.get('cpus')
    threads_arg = "" if threads is None else f"-t {threads}"
    # matches megahit's convention: pin to the allocation we actually got
    # (85% of it), not spades' own node-wide default.
    mem_gb = context.params.get('memory')
    mem_arg = f"-m {max(1, int(mem_gb * 0.85))}" if mem_gb else ""

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            spades.py {mode} {threads_arg} {mem_arg} \
                {reads_arg} \
                -o spades_ws
            [[ $(head spades_ws/contigs.fasta | wc --chars) -ne 0 ]] && mv spades_ws/contigs.fasta {iout.container} || echo "assembly was empty"
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
        cpus=8,
        memory=Size.GB(32),
        duration=Duration(hours=18),
    )
)
