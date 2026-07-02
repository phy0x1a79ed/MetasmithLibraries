"""ORF-call the reference fosmid inserts (prodigal) -> ORF proteins + GFF.

Fresh ORF calling on the 132 >=29kb reference inserts, source-tagged as the
`fosmid` ORFs the four ECSPr evidence lanes annotate. Mamba-native (prodigal in
the fabfos-annot env) rather than the docker pprodigal wrapper in metagenomics/,
so it runs under Runtime.MAMBA. `-p meta` (metagenomic mode) since each insert is
a short standalone contig. The amino-acid ORF FASTA IS sequences::open_reading_frames,
the query for kofam / dl_ec / uniref / embed lanes.
"""
from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
exp     = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
inserts = model.AddRequirement(lib.GetType("fosmids::reference_inserts"), parents={exp})
env     = model.AddRequirement(lib.GetType("envs::prodigal.condaenv"))
cds     = model.AddProduct(lib.GetType("sequences::open_reading_frames"))
gff     = model.AddProduct(lib.GetType("sequences::gff"))

def protocol(context: ExecutionContext):
    iins = context.Input(inserts)
    icds = context.Output(cds)
    igff = context.Output(gff)
    context.ExecWithContainer(
        image=env,
        cmd=f"""prodigal \
            -p meta \
            -i {iins.container} \
            -a {icds.container} \
            -f gff \
            -o {igff.container}""",
    )
    return ExecutionResult(
        manifest=[{cds: icds.local, gff: igff.local}],
        success=icds.local.exists() and igff.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=2, memory=Size.GB(4), duration=Duration(hours=1)),
)
