"""CheckM2 (drop-in replacement for the original CheckM v1 transform).

Public surface — model name, dtype produced (taxonomy::checkm_stats), and the
per-genome TSV shape consumed by the aggregator — is unchanged. Internals
swapped to `checkm2 predict`; the CheckM2 DiamondDB is baked into the image
and exposed via CHECKM2DB, so no --database_path is needed.

The aggregator (_parse_checkm) joins on a "Bin Id" column. CheckM2's
quality_report.tsv calls that column "Name"; we rename it on emit.
"""
import csv
import shutil
from pathlib import Path
from metasmith.python_api import *

lib         = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model       = Transform()
image       = model.AddRequirement(lib.GetType("env::checkm.env")) # checkm2; DiamondDB baked, CHECKM2DB set
asm         = model.AddRequirement(lib.GetType("sequences::putative_genome"))
out         = model.AddProduct(lib.GetType("taxonomy::checkm_stats"))


def protocol(context: ExecutionContext):
    input_dir = Path("input")
    input_dir.mkdir()

    in2out = {}  # asm_stem -> iout
    for item in context.AsBatch():
        iasm = item.Input(asm)
        iout = item.Output(out)
        in2out[iasm.local.stem] = iout
        src = iasm.local
        dest = input_dir / iasm.local.name
        Log.Info(f"registering genome [{src}] -> [{dest}]")
        shutil.copy(src, dest, follow_symlinks=True)

    threads = context.params.get('cpus')
    threads = "" if threads is None else f"--threads {threads}"
    ext = iasm.container.suffix.replace(".", "")
    out_dir = "checkm2_out"
    report = f"{out_dir}/quality_report.tsv"
    context.ExecWithContainer(
        image = image,
        cmd = f"""
            export PATH=/opt/conda/envs/external_checkm2_env/bin:/opt/conda/bin:$PATH
            checkm2 predict {threads} -x {ext} --input ./input --output-directory ./{out_dir}
        """
    )

    # Split combined quality_report.tsv → one per-genome TSV.
    # Rename CheckM2's "Name" column to "Bin Id" so the aggregator's
    # _parse_checkm() (which keys on Bin Id) keeps working.
    with open(report, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        name_col = header.index("Name")
        out_header = header.copy()
        out_header[name_col] = "Bin Id"
        rows_by_stem = {}
        for row in reader:
            if not row:
                continue
            rows_by_stem[row[name_col]] = row

    manifest = []
    for asm_stem, iout in in2out.items():
        with open(iout.local, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(out_header)
            row = rows_by_stem.get(asm_stem)
            if row is not None:
                w.writerow(row)
        manifest.append({out: iout.local})

    return ExecutionResult(
        manifest=manifest,
        success=Path(report).exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    batch_size=25,
    resources=Resources(
        cpus=2,
        memory=Size.GB(16),  # checkm2 is lighter than v1 (~12 GB peak); revisit after smoke
        duration=Duration(hours=1),
    )
)
