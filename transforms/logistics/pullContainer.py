from pathlib import Path
from metasmith.python_api import *

lib     = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model   = Transform()
image   = model.AddRequirement(lib.GetType("containers::container"))
log     = model.AddProduct(lib.GetType("containers::pulled_container"))

def protocol(context: ExecutionContext):
    ilog=context.Output(log)
    container = context.GetContainerModel(image)
    success = Path("pull_success")
    context.external_shell.Exec(f"[[ -e {container.GetLocalPath()} ]] && touch {success}")
    if success.exists():
        with open(ilog.local, "w") as f:
            f.write(f"local image already exists at [{container.GetLocalPath()}]")
    else:
        context.external_shell.Exec(f'bash -c "{container.MakePullCommand()} 2>&1 && touch {success}" | tee {ilog.external}')
    return ExecutionResult(
        manifest=[{log: ilog.local}],
        success=success.exists()
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(
        cpus=1,
        memory=Size.GB(1),
        duration=Duration(hours=1),
    )
)
