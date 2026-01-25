import os, sys
from pathlib import Path
from metasmith.python_api import Agent, Source, SshSource, DataInstanceLibrary, TransformInstanceLibrary, DataTypeLibrary
from metasmith.python_api import Resources, Size, Duration
from metasmith.python_api import ContainerRuntime

base_dir = Path("./cache")

# agent_home = Source.FromLocal((base_dir/"local_home").resolve())
agent_home = Source.FromLocal("/home/tony/workspace/tools/Metasmith/main/local_mock/cache/local_home")
smith = Agent(
    home = agent_home,
    # runtime=ContainerRuntime.APPTAINER,
    runtime=ContainerRuntime.DOCKER,
)

# agent_home = SshSource(host="chamois", path=Path("/home/tliu/metasmith")).AsSource()
# smith = Agent(
#     home = agent_home,
#     runtime=ContainerRuntime.APPTAINER,
#     setup_commands=[
#         'PATH=/home/tliu/miniforge3/envs/def/bin/:$PATH',
#         'export TMPDIR="/home/$USER/tmp"',
#         'mkdir -p $TMPDIR',
#         'export APPTAINER_CACHEDIR="$TMPDIR"',
#         'export APPTAINER_TMPDIR="$TMPDIR"',
#     ]
# )

# agent_home = SshSource(host="sockeye", path=Path("/scratch/st-shallam-1/pwy_group/metasmith")).AsSource()
# smith = Agent(
#     home = agent_home,
#     runtime=ContainerRuntime.APPTAINER,
#     setup_commands=[
#         'module load gcc/9.4.0',
#         'module load apptainer/1.3.1',
#         'export APPTAINER_CACHEDIR="/scratch/st-shallam-1/pwy_group/apptainer_cache"',
#     ]
# )

# agent_home = SshSource(host="fir", path=Path("/scratch/phyberos/metasmith")).AsSource()
# smith = Agent(
#     home = agent_home,
#     runtime=ContainerRuntime.APPTAINER,
#     setup_commands=[
#         "module load StdEnv/2023",
#         "module load apptainer/1.3.5",
#     ]
# )

# smith.Deploy()

# import ipynbname
# notebook_name = ipynbname.name()
notebook_name = Path(__file__).stem
in_dir = base_dir/f"{notebook_name}/inputs.xgdb"

inputs = DataInstanceLibrary(in_dir)
inputs.Purge()
inputs.AddTypeLibrary("ncbi", DataTypeLibrary.Load("../data_types/ncbi.yml"))
inputs.AddTypeLibrary("sequences", DataTypeLibrary.Load("../data_types/sequences.yml"))
inputs.AddTypeLibrary("pangenome", DataTypeLibrary.Load("../data_types/pangenome.yml"))

group = inputs.AddValue("pangenome", "e coli", "pangenome::pangenome")
inputs.AddValue("DH10b", "GCF_000019425.1", "ncbi::assembly_accession", parents={group})
inputs.AddValue("K12", "GCF_000005845.2", "ncbi::assembly_accession", parents={group})
inputs.AddItem((base_dir/f"{notebook_name}/epi300.gbk").resolve(), "sequences::gbk", parents={group})
inputs.LocalizeContents()
inputs.Save()

# inputs = DataInstanceLibrary.Load(in_dir)

resources = [
    DataInstanceLibrary.Load(f"../resources/{n}")
    for n in ["containers", "lib"]
]

transforms = [
    TransformInstanceLibrary.Load(f"../transforms/{n}")
    for n in ["logistics", "pangenome"]
]

task = smith.GenerateWorkflow(
    samples=inputs.AsSamples(),
    resources=resources,
    transforms=transforms,
    # targets=[inputs.GetType("sequences::gbk")]
    targets=["pangenome::heatmap"]
)
task.plan.RenderDAG(base_dir/f"{notebook_name}/dag")
print(task.ok, len(task.plan.steps))

smith.StageWorkflow(task, on_exist="clear")
params = dict(
    executor=dict(
        cpus=14,
    ),
    process=dict(
        tries=1,
    ),
)
smith.RunWorkflow(
    task=task,
    config_file=smith.GetNxfConfigPresets()["local"],
    params=params,
    resource_overrides={
        "*": Resources(
            memory=Size.GB(1),
            cpus=14,
        )
    }
)

# # with open("../secrets/slurm_account_fir") as f:
# #     SLURM_ACCOUNT = f.readline()
# # smith.RunWorkflow(
# #     task,
# #     config_file=smith.GetNxfConfigPresets()["slurm"],
# #     params=dict(
# #         slurmAccount=SLURM_ACCOUNT,
# #     )
# # )