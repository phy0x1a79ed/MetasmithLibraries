from metasmith.python_api import *
from pathlib import Path

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image_pbert = model.AddRequirement(lib.GetType("containers::proteinbert.oci"))
image_polars = model.AddRequirement(lib.GetType("containers::polars.oci"))
orfs = model.AddRequirement(lib.GetType("sequences::orfs"))
out_embeddings = model.AddProduct(lib.GetType("annotation::proteinbert_embeddings"))
out_index = model.AddProduct(lib.GetType("annotation::proteinbert_index"))


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    iemb = context.Output(out_embeddings)
    iidx = context.Output(out_index)

    threads = context.params.get("cpus", 8)

    # Run ProteinBERT to generate embeddings
    context.ExecWithContainer(
        image=image_pbert,
        cmd=f"""
            pbert run \
                -i {iorfs.container} \
                -o pbert_output \
                --threads {threads} \
                --protein_size 512 \
                --model_batch 2048 \
                -x 1
        """,
    )

    # Create combiner script using polars
    combiner_script = Path("combine_embeddings.py")
    with open(combiner_script, "w") as f:
        f.write('''
import sys
import numpy as np
import polars as pl
from pathlib import Path

input_dir = Path(sys.argv[1])
output_file = Path(sys.argv[2])

npy_files = sorted(input_dir.glob("*.npy"))
if npy_files:
    arrays = [np.load(f) for f in npy_files]
    combined = np.vstack(arrays)[:, -512:]
    col_names = [f"dim_{i}" for i in range(512)]
    df = pl.DataFrame(combined, schema=col_names)
    df.write_parquet(output_file)
''')

    # Combine embeddings using polars container
    context.ExecWithContainer(
        image=image_polars,
        cmd=f"""
            python {combiner_script} pbert_output {iemb.container}
        """,
    )

    # Copy index file
    index_files = list(Path("pbert_output").glob("*.csv"))
    if index_files:
        context.LocalShell(f"cp {index_files[0]} {iidx.local}")
    else:
        with open(iidx.local, "w") as f:
            f.write("sequence_id,index\n")

    return ExecutionResult(
        manifest=[
            {
                out_embeddings: iemb.local,
                out_index: iidx.local,
            },
        ],
        success=iemb.local.exists() and iidx.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=8,
        memory=Size.GB(32),
        duration=Duration(hours=12),
    ),
)
