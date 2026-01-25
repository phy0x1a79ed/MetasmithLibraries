import os, sys
from pathlib import Path
from metasmith.python_api import Agent, Source, SshSource, DataInstanceLibrary, TransformInstanceLibrary, DataTypeLibrary
from metasmith.python_api import Resources, Size, Duration
from metasmith.python_api import ContainerRuntime

MLIB = Path("/home/tony/workspace/tools/MetasmithLibraries")

base_dir = Path("./cache")

# agent_home = Source.FromLocal((base_dir/"local_home").resolve())
# smith = Agent(
#     home = agent_home,
#     runtime=ContainerRuntime.APPTAINER,
#     # runtime=ContainerRuntime.DOCKER,
# )

agent_home = SshSource(host="sockeye", path=Path("/scratch/st-shallam-1/pwy_group/metasmith")).AsSource()
smith = Agent(
    home = agent_home,
    runtime=ContainerRuntime.APPTAINER,
    setup_commands=[
        'module load gcc/9.4.0',
        'module load apptainer/1.3.1',
    ]
)

# agent_home = SshSource(host="fir", path=Path("/scratch/phyberos/metasmith")).AsSource()
# smith = Agent(
#     home = agent_home,
#     runtime=ContainerRuntime.APPTAINER,
#     setup_commands=[
#         'module load gcc/9.4.0',
#         'module load apptainer/1.3.1',
#     ]
# )

# smith.Deploy(assertive=True)
# sys.exit(0)

notebook_name = Path(__file__).stem
in_dir = base_dir/f"{notebook_name}/inputs.xgdb"

containers = DataInstanceLibrary.Load(MLIB/"resources/containers")
logistics = TransformInstanceLibrary.Load(MLIB/f"transforms/logistics")

task = smith.GenerateWorkflow(
    samples=containers.AsSamples(),
    resources=[],
    transforms=[logistics],
    # targets=[inputs.GetType("sequences::gbk")]
    targets=["containers::pulled_container"]
)
# task.plan.RenderDAG("pull_dag.svg", blacklist_namespaces=set())
print(task.ok, len(task.plan.steps))

# smith.StageWorkflow(task, on_exist="clear")
smith.StageWorkflow(task, on_exist="update_all")
# smith.StageWorkflow(task, on_exist="update_workflow")

params = dict(
    executor=dict(
        cpus=15,
        memory='20 GB',
        queueSize=5,
    ),
)
smith.RunWorkflow(
    task,
    config_file=smith.GetNxfConfigPresets()["local"],
    params=params,
    # resource_overrides={
    #     "all": Resources(
    #         memory=Size.GB(2),
    #         cpus=Size.GB(2),
    #     )   
    # }
)
