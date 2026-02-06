from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
src     = model.AddRequirement(lib.GetType("functional_annotation::protein_ref_url"))
fasta   = model.AddProduct(lib.GetType("functional_annotation::protein_ref_fasta"))

def protocol(context: ExecutionContext):
    isrc   = context.Input(src)
    ifasta = context.Output(fasta)

    with open(isrc.local) as f:
        url = f.readline().strip()

    is_gzipped = url.endswith(".gz")
    dl_name = "ref_download.gz" if is_gzipped else "ref_download.faa"

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            python3 -c "
import urllib.request, sys
print('Downloading protein reference from {url} ...')
urllib.request.urlretrieve('{url}', '{dl_name}')
print('Done.')
"
            {'gunzip -c ' + dl_name + ' > ' + str(ifasta.container) if is_gzipped else 'mv ' + dl_name + ' ' + str(ifasta.container)}
        """,
    )

    return ExecutionResult(
        manifest=[
            {
                fasta: ifasta.local,
            },
        ],
        success=ifasta.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=src,
    labels=["local"],
    resources=Resources(
        cpus=1,
        memory=Size.GB(2),
        duration=Duration(hours=4),
    )
)
