from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
img_ipr = model.AddRequirement(lib.GetType("containers::interproscan.oci"))
data    = model.AddProduct(lib.GetType("annotation::interproscan_data"))

# InterProScan 5.67-99.0 data bundle
IPRSCAN_DATA_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz"

def protocol(context: ExecutionContext):
    idata = context.Output(data)

    context.ExecWithContainer(
        image=image,
        cmd=f"""
            wget -q {IPRSCAN_DATA_URL} -O interproscan-data.tar.gz
            mkdir -p {idata.container}
            tar xzf interproscan-data.tar.gz -C {idata.container} --strip-components=1
        """,
    )

    # compile hmmer indexes
    context.ExecWithContainer(
        image=img_ipr,
        binds=[(idata.external/"data", "/opt/interproscan/data")],
        cmd=f"""\
            python3 setup.py -f interproscan.properties --force
        """
    )

    return ExecutionResult(
        manifest=[{data: idata.local}],
        success=idata.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=img_ipr,
    labels=["local"],
    resources=Resources(
        cpus=1,
        memory=Size.GB(8),
        duration=Duration(hours=8),
    ),
)
