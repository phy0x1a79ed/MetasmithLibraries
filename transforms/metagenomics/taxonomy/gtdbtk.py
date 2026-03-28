from metasmith.python_api import *
from pathlib import Path
import shutil

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("containers::gtdbtk.oci"))
ref         = model.AddRequirement(lib.GetType("taxonomy::gtdb"))
asm         = model.AddRequirement(lib.GetType("sequences::assembly"))
tax         = model.AddProduct(lib.GetType("taxonomy::gtdbtk"))
# raw         = model.AddProduct(lib.GetType("taxonomy::gtdbtk_raw"))

def protocol(context: ExecutionContext):
    iref    = context.Input(ref)
    # iraw    = context.Output(raw)

    genome_dir = Path("./assemblies")
    genome_dir.mkdir()
    in2out = {}
    for item in context.AsBatch():
        iasm    = item.Input(asm)
        itax    = item.Output(tax)
        in2out[iasm.local.stem] = itax.local.name # gtdb removes the file extension, so .stem
        src = iasm.local
        dest = genome_dir/iasm.local.name
        Log.Info(f"registering genome [{src}] -> [{dest}]")
        shutil.copy(src, dest, follow_symlinks=True)
    
    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--cpus {threads}"
    mem = context.params.get('memory')
    if mem:
        _mem_gb = int(float(mem))
        pplacer_cpus = f"--pplacer_cpus {max(1, (_mem_gb-8)//40)}"
    else:
        pplacer_cpus = ""

    ext = iasm.container.suffix.replace(".", "")
    TEMP_PREFIX = "temp"
    # temp_scratch = Path(f"{TEMP_PREFIX}.scratch") # replaces RAM
    temp_ws = Path(f"{TEMP_PREFIX}.ws")
    # reduce pplacer memory usage by writing to disk (slower).
    # --scratch_dir {temp_scratch} \
    # Skip the skani ANI screening step to classify genomes.
    # --skip_ani_screen
    # todo:
    # - use prodigal genes
    # - separate the skani screen? is this needed for small runs? 
    # - batchify
    out_raw = Path("./gtdb_raw")
    context.ExecWithContainer(
        binds=[
            (iref.external, "/ref"),
        ],
        image = image,
        cmd = f"""\
            mkdir -p {temp_ws}
            export GTDBTK_DATA_PATH=/ref
            gtdbtk classify_wf -x {ext} \
                {threads} {pplacer_cpus} \
                --force \
                --skip_ani_screen \
                --tmpdir {temp_ws} \
                --genome_dir {genome_dir} \
                --out_dir {out_raw}
        """
    )
    
    file_candidates = [p for p in out_raw.glob("classify/*summary.tsv")]
    assert len(file_candidates)>0, "no *summary.tsv files"
    rows = {}
    for table in file_candidates:
        with open(table) as tsv:
            header = tsv.readline()
            for l in tsv:
                toks = l.strip().split("\t")
                k = toks[0]
                rows[k] = l, header

    # have output file creation order match input order, since paranoid
    for _asm, _tax in in2out.items():
        row, header = rows[_asm]
        if row[:-1] != "\n": row+="\n"
        if header[:-1] != "\n": header+="\n"
        with open(_tax, "w") as out:
            out.write(header)
            out.write(row)

    # todo: return trees as well
    # out_tree = context.output_folder.joinpath(f"{sample}.tax.tree")
    # file_candidates = os.listdir(classify_out) if classify_out.exists() else []
    # file_candidates = [classify_out.joinpath(f) for f in file_candidates if (f.endswith(".tree") and "backbone" not in str(f))]
    # with open(out_tree, "w") as out:

    return ExecutionResult(
        manifest=[
            {
                tax: itax.local,
            },
        ],
        success=itax.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    batch_size=100,
    resources=Resources(
        cpus=2,
        memory=Size.GB(120), # r226 used 107 GB
        duration=Duration(hours=12),
    )
)
