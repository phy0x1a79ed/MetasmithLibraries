from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
src     = model.AddRequirement(lib.GetType("functional_annotation::kofam_source"))
db      = model.AddProduct(lib.GetType("functional_annotation::kofam_db"))

KOFAM_PROFILES_URL = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
KOFAM_KO_LIST_URL  = "https://www.genome.jp/ftp/db/kofam/ko_list.gz"

def protocol(context: ExecutionContext):
    idb = context.Output(db)

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            mkdir -p {idb.container}
            python3 -c "
import urllib.request, sys
print('Downloading KOfam profiles...')
urllib.request.urlretrieve('{KOFAM_PROFILES_URL}', 'profiles.tar.gz')
print('Downloading KOfam ko_list...')
urllib.request.urlretrieve('{KOFAM_KO_LIST_URL}', 'ko_list.gz')
print('Done.')
"
            tar xzf profiles.tar.gz -C {idb.container}
            gunzip -c ko_list.gz > {idb.container}/ko_list
        """,
    )

    return ExecutionResult(
        manifest=[
            {
                db: idb.local,
            },
        ],
        success=(idb.local / "profiles").exists() and (idb.local / "ko_list").exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=src,
    labels=["local"],
    resources=Resources(
        cpus=1,
        memory=Size.GB(2),
        duration=Duration(hours=2),
    )
)
