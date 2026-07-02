"""Lane 1 -- kofamscan KO assignments for the fosmid ORFs -> kofam_hits.

Fresh kofamscan run (mamba, fabfos-kofam env) over the reference-insert ORFs
against the staged KOfam HMM profiles + ko_list. Emits an above-threshold CSV in
the column shape `resources/lib/fabfos_evidence.py::read_kofam` consumes
(orf, ko, hmm_threshold, score, description); the compiler projects KO -> MNXR.

GATED: needs the KOfam HMM DB (profiles dir ~4GB + ko_list) staged and the
fabfos-kofam conda env (kofamscan + hmmer) provisioned.
"""
from pathlib import Path
from metasmith.python_api import *

lib      = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model    = Transform()
exp      = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
orfs     = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
profiles = model.AddRequirement(lib.GetType("functional_annotation::kofam_profiles"))
ko_list  = model.AddRequirement(lib.GetType("functional_annotation::kofam_ko_list"))
env      = model.AddRequirement(lib.GetType("envs::kofamscan.condaenv"))
hits     = model.AddProduct(lib.GetType("functional_annotation::kofam_hits"))

def _parse(detail_txt, out_csv):
    # kofamscan --format=detail: '*'? gene KO thrshld score E-value "KO definition"
    with open(detail_txt) as fin, open(out_csv, "w") as fout:
        fout.write("orf,ko,hmm_threshold,score,description\n")
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.lstrip("* ").rstrip("\n").split(None, 5)
            if len(parts) < 5:
                continue
            gene, ko, thr, score, _ev = parts[:5]
            desc = parts[5].replace(",", ";") if len(parts) > 5 else ""
            try:
                if float(score) >= float(thr):
                    fout.write(f"{gene},{ko},{thr},{score},{desc}\n")
            except ValueError:
                continue

def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    iprof = context.Input(profiles)
    iko   = context.Input(ko_list)
    ohits = context.Output(hits)
    cpus  = context.params.get("cpus", 8)
    context.ExecWithContainer(
        image=env,
        cmd=f"""exec_annotation \
            -o kofam_detail.txt \
            --profile={iprof.container} \
            --ko-list={iko.container} \
            --cpu={cpus} \
            --e-value=0.01 \
            --format=detail \
            --no-report-unannotated \
            {iorfs.container}""",
    )
    _parse("kofam_detail.txt", str(ohits.local))
    return ExecutionResult(
        manifest=[{hits: ohits.local}],
        success=ohits.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=8, memory=Size.GB(16), duration=Duration(hours=8)),
)
