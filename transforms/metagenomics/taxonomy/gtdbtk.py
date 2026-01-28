from metasmith.python_api import *
from pathlib import Path

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("containers::gtdbtk.oci"))
ref         = model.AddRequirement(lib.GetType("taxonomy::gtdb"))
asm         = model.AddRequirement(lib.GetType("sequences::assembly"))
tax         = model.AddProduct(lib.GetType("taxonomy::gtdbtk"))
raw         = model.AddProduct(lib.GetType("taxonomy::gtdbtk_raw"))

def protocol(context: ExecutionContext):
    iref    = context.Input(ref)
    iasm    = context.Input(asm)
    itax    = context.Output(tax)
    iraw    = context.Output(raw)
    
    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--cpus {threads}"
    mem = context.params.get('memory')
    if mem:
        _mem_gb = int(float(mem))
        pplacer_cpus = f"--pplacer_cpus {max(1, (_mem_gb-4)//40)}"
    else:
        pplacer_cpus = ""

    ext = iasm.container.suffix.replace(".", "")
    TEMP_PREFIX = "temp"
    temp_scratch = Path(f"{TEMP_PREFIX}.scratch")
    temp_ws = Path(f"{TEMP_PREFIX}.ws")
    # reduce pplacer memory usage by writing to disk (slower).
    # --scratch_dir {temp_scratch} \
    # Skip the skani ANI screening step to classify genomes.
    # --skip_ani_screen
    # todo:
    # - use prodigal genes
    # - separate the skani screen? is this needed for small runs? 
    # - batchify
    context.ExecWithContainer(
        binds=[
            (iref.external, "/ref"),
        ],
        image = image,
        cmd = f"""\
            mkdir -p input {temp_scratch} {temp_ws}
            cp -L {iasm.container} ./input/
            export GTDBTK_DATA_PATH=/ref
            gtdbtk classify_wf -x {ext} \
                {threads} {pplacer_cpus} \
                --force \
                --skip_ani_screen \
                --tmpdir {temp_ws} \
                --genome_dir ./input \
                --out_dir {iraw.container}
        """
    )
    
    file_candidates = [p for p in iraw.local.glob("classify/*summary.tsv")]
    assert len(file_candidates)>0, "no *summary.tsv files"
    _written_header = False
    with open(itax.local, "w") as out:
        for table in file_candidates:
            with open(table) as tsv:
                header = tsv.readline()
                if not _written_header: out.write(header); _written_header = True
                for l in tsv:
                    out.write(l)

    # todo: return trees as well
    # out_tree = context.output_folder.joinpath(f"{sample}.tax.tree")
    # file_candidates = os.listdir(classify_out) if classify_out.exists() else []
    # file_candidates = [classify_out.joinpath(f) for f in file_candidates if (f.endswith(".tree") and "backbone" not in str(f))]
    # with open(out_tree, "w") as out:

    return ExecutionResult(
        manifest=[
            {
                tax: itax.local,
                raw: iraw.local,
            },
        ],
        success=itax.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=1,
        memory=Size.GB(80),
        duration=Duration(hours=4),
    )
)
