from metasmith.python_api import *
from pathlib import Path
import csv

lib = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("containers::interproscan.oci"))
orfs = model.AddRequirement(lib.GetType("sequences::orfs"))
data_dir = model.AddRequirement(lib.GetType("annotation::interproscan_data"))
out_gff = model.AddProduct(lib.GetType("annotation::interproscan_results"))


def parse_interpro_gff(input_path, output_path):
    headers = ["orf", "start", "stop", "score", "database", "Name", "Dbxref"]
    with open(input_path, 'r') as gff_file, open(output_path, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != "protein_match":
                continue
            score = parts[5]
            score = "" if score == "." else score
            attrs = {}
            for item in parts[8].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    attrs[key] = value.strip('"')
            dbxref = attrs.get("Dbxref", "")
            if dbxref.startswith("InterPro:"):
                dbxref = dbxref.replace("InterPro:", "", 1)
            writer.writerow({
                "orf": parts[0],
                "start": parts[3],
                "stop": parts[4],
                "score": score,
                "database": parts[1],
                "Name": attrs.get("Name", ""),
                "Dbxref": dbxref,
            })


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    idata = context.Input(data_dir)
    ijson = context.Output(out_json)
    igff = context.Output(out_gff)

    cpus = context.params.get("cpus")
    cpus_string = "" if cpus is None else f"--cpu {cpus}"

    # sed to remove stop codon asterisks from protein sequences
    # Run InterProScan with external data directory
    context.ExecWithContainer(
        image=image,
        binds=[(idata.external/"data", "/opt/interproscan/data")],
        cmd=f"""
            mkdir -p output
            sed '/^[^>]/s/\*//g' {iorfs.container} > ./{iorfs.container.stem}.clean.faa
            /opt/interproscan/interproscan.sh \
                --disable-precalc \
                --verbose \
                --seqtype p \
                {cpus_string} \
                -i ./{iorfs.container.stem}.clean.faa \
                -f gff3 \
                -d output
        """,
    )

    # Copy JSON output
    json_files = list(Path("output").glob("*.json"))
    if json_files:
        context.LocalShell(f"cp {json_files[0]} {ijson.local}")

    # Parse GFF3 output into CSV
    gff_files = list(Path("output").glob("*.gff3"))
    if gff_files:
        parse_interpro_gff(str(gff_files[0]), str(igff.local))

    return ExecutionResult(
        manifest=[
            {
                out_json: ijson.local,
                out_gff: igff.local,
            },
        ],
        success=ijson.local.exists() and igff.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=16,
        memory=Size.GB(32),
        duration=Duration(hours=24),
    ),
)
