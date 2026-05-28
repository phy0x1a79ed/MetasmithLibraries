from metasmith.python_api import *
from pathlib import Path

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("containers::promotech.oci"))
assembly   = model.AddRequirement(lib.GetType("sequences::assembly"))
gff        = model.AddRequirement(lib.GetType("sequences::gff"))
extractor  = model.AddRequirement(lib.GetType("lib::extract_noncoding_chunks.py"))
merger     = model.AddRequirement(lib.GetType("lib::merge_promotech_results.py"))

out_pred = model.AddProduct(lib.GetType("annotation::promotech_predictions"))


def _resolve(path):
    """Return container path, handling both absolute and relative cases."""
    return path.container if path.container.is_absolute() else Path("/ws") / path.container


def protocol(context: ExecutionContext):
    iasm = context.Input(assembly)
    igff = context.Input(gff)
    iextract = context.Input(extractor)
    imerge = context.Input(merger)
    opred = context.Output(out_pred)

    # /ws is the writable workspace inside the container (Nextflow work dir)
    chunks_dir = "/ws/pt_chunks"
    results_dir = "/ws/pt_results"

    # Step 1: Extract non-coding intervals into ≤1Mbp chunk FASTAs
    context.ExecWithContainer(
        image=image,
        cmd=f"""
            python {_resolve(iextract)} \
                --fasta {_resolve(iasm)} \
                --gff {_resolve(igff)} \
                --outdir {chunks_dir} \
                --max-chunk-bp 1000000 &&
            cd /ws
        """,
    )

    # Step 2: Run PromoTech parse + predict on each chunk in parallel
    # xargs -P runs N chunks concurrently (1 thread each)
    cpus = context.params.get("cpus", 4)
    context.ExecWithContainer(
        image=image,
        cmd=f"""
            cd /opt/promotech &&
            mkdir -p {results_dir} &&
            run_chunk() {{
                chunk=$1
                name=$(basename "$chunk" .fna)
                mkdir -p {results_dir}/$name
                python promotech.py -pg -f "$chunk" \
                    -o {results_dir}/$name -m RF-HOT 2>&1 &&
                python promotech.py -g -i {results_dir}/$name \
                    -o {results_dir}/$name -m RF-HOT -t 0.5 2>&1
            }} &&
            export -f run_chunk &&
            ls {chunks_dir}/*.fna | xargs -I{{}} -P {cpus} \
                bash -c 'run_chunk "$@"' _ {{}} &&
            cd /ws
        """,
    )

    # Step 3: Merge results and remap coordinates
    context.ExecWithContainer(
        image=image,
        cmd=f"""
            python {_resolve(imerge)} \
                --manifest {chunks_dir}/manifest.json \
                --results-dir {results_dir} \
                --output /ws/pt_merged.tsv &&
            cd /ws
        """,
    )

    context.LocalShell(f"cp /ws/pt_merged.tsv {opred.local}")

    return ExecutionResult(
        manifest=[{out_pred: opred.local}],
        success=opred.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=assembly,
    resources=Resources(
        cpus=4,
        memory=Size.GB(64),
        duration=Duration(hours=4),
    ),
)
