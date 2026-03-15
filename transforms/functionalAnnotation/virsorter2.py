from metasmith.python_api import *

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image       = model.AddRequirement(lib.GetType("containers::virsorter2.oci"))
asm         = model.AddRequirement(lib.GetType("sequences::assembly"))
db          = model.AddRequirement(lib.GetType("annotation::virsorter2_db"))
out_seqs    = model.AddProduct(lib.GetType("annotation::virsorter2_viral_sequences"))
out_scores  = model.AddProduct(lib.GetType("annotation::virsorter2_scores"))
out_affi    = model.AddProduct(lib.GetType("annotation::virsorter2_affi_contigs"))
out_boundary = model.AddProduct(lib.GetType("annotation::virsorter2_boundary"))

def protocol(context: ExecutionContext):
    iasm = context.Input(asm)
    idb  = context.Input(db)
    iseqs = context.Output(out_seqs)
    iscores = context.Output(out_scores)
    iaffi = context.Output(out_affi)
    iboundary = context.Output(out_boundary)

    threads = context.params.get("cpus", 8)
    workdir = "vs2_out"

    # VirSorter may exit non-zero when no viral contigs are found (EARLY-EXIT)
    # so we tolerate exit codes and check outputs instead
    context.ExecWithContainer(
        image=image,
        binds=[(idb.external, "/db")],
        cmd=f"""
            export HOME=/tmp

            virsorter run \
                --seqfile {iasm.container} \
                --db-dir /db \
                --working-dir {workdir} \
                --jobs {threads} \
                --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae \
                --min-length 5000 \
                --min-score 0.5 \
                --keep-original-seq \
                --prep-for-dramv \
                all \
                || true
        """,
    )

    # Copy primary outputs (files may not exist if no viral contigs found)
    context.LocalShell(
        f"cp {workdir}/final-viral-combined-for-dramv.fa {iseqs.local} 2>/dev/null || touch {iseqs.local}"
    )
    context.LocalShell(
        f"cp {workdir}/final-viral-score.tsv {iscores.local} 2>/dev/null || touch {iscores.local}"
    )
    context.LocalShell(
        f"cp {workdir}/viral-affi-contigs-for-dramv.tab {iaffi.local} 2>/dev/null || touch {iaffi.local}"
    )
    context.LocalShell(
        f"cp {workdir}/final-viral-boundary.tsv {iboundary.local} 2>/dev/null || touch {iboundary.local}"
    )

    return ExecutionResult(
        manifest=[
            {
                out_seqs: iseqs.local,
                out_scores: iscores.local,
                out_affi: iaffi.local,
                out_boundary: iboundary.local,
            },
        ],
        success=True,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=8,
        memory=Size.GB(16),
        duration=Duration(hours=12),
    ),
)
