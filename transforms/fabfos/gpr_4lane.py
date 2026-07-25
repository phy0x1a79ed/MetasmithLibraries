"""GPR mapper (chosen 4 lanes) -- fold kofam + CLEAN + uniref50 + ProteinBERT into
one gene-attributed GPR table -> annotation::gpr_table_4lane.

The "native run tool -> mapper into GPR table format" half of the pipeline: it
takes the four lane outputs and translates each lane's identifiers to MetaNetX
reactions (KO->MNXR, EC->MNXR, UniProt->MNXR) and does kNN label transfer for the
ProteinBERT embedding lane, emitting the frozen long-format 8-col GPR schema:

    source, orf, channel, mnxr, intermediate_id, intermediate_name, raw_score, projection_via

one row per (ORF, lane channel, reaction). No cross-lane dedup -- the same
ORF->MNXR claim from two lanes stays two rows, distinguished by `channel`.

Modeled on ptools_annotation_gather.py (require orfs + the lane outputs,
group_by=orfs, python_for_data_science image, long-format parquet). The bridge
translation reuses resources/lib/fabfos_evidence.py (bridge loaders + stitle
cleaner + schema + score floors, pandas-only); the embedding kNN vote is a
lightweight numpy port of resources/lib/fabfos_embed_transfer.py::apply_one.

Chosen-4 = one lane per identifier space plus the current best EC lane:
  kofam (KO->MNXR), CLEAN (EC->MNXR), uniref50 (UniProt->MNXR via BSR),
  ProteinBERT (embedding kNN transfer).

GATED: the KO/EC/UniProt->MNXR bridges and the reference_label_pool are staged
references (their build transforms are deferred); planning is type-driven, so this
resolves with none of them present. Running it end-to-end is a later milestone.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
orfs      = model.AddRequirement(lib.GetType("sequences::orfs"))
kofam     = model.AddRequirement(lib.GetType("annotation::kofamscan_results"))
clean     = model.AddRequirement(lib.GetType("annotation::clean_predictions"))
uniref    = model.AddRequirement(lib.GetType("annotation::diamond_uniref50_results"))
pbert_emb = model.AddRequirement(lib.GetType("annotation::proteinbert_embeddings"))
pbert_idx = model.AddRequirement(lib.GetType("annotation::proteinbert_index"))
ko_br     = model.AddRequirement(lib.GetType("ref::ko_to_mnxr"))
ec_br     = model.AddRequirement(lib.GetType("ref::ec_to_mnxr"))
up_br     = model.AddRequirement(lib.GetType("ref::uniprot_to_mnxr"))
pool      = model.AddRequirement(lib.GetType("ref::reference_label_pool"))
ev_lib    = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
knn_lib   = model.AddRequirement(lib.GetType("lib::fabfos_embed_transfer.py"))
out_gpr   = model.AddProduct(lib.GetType("annotation::gpr_table_4lane"))


# The inline driver reads the dev2-format lane outputs, projects identifiers to
# MNXR via the staged bridges, runs the embedding kNN transfer against the
# reference pool, and writes the 8-col GPR parquet. `{...}` placeholders are
# filled from the container paths; `{{`/`}}` are literal braces.
DRIVER = r'''
import sys, os
import numpy as np
import pandas as pd

# reuse the ported bridge loaders / stitle cleaner / schema / floors (pandas-only)
sys.path.insert(0, os.path.dirname("{ev_lib}"))
import fabfos_evidence as fe

SCHEMA = fe.SCHEMA_COLS
SOURCE = "fosmid"
K = 30
PBERT_FLOOR = 0.20

def orf_ids(fasta):
    out = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].split()[0])
    return out

# ---- kofam lane: dev2 kofamscan_results = gene_name,KO,thrshld,score,E-value,best
def lane_kofam(path, ko_to_mnxr):
    df = pd.read_csv(path)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["thrshld"] = pd.to_numeric(df["thrshld"], errors="coerce")
    df = df[df["score"].notna() & df["thrshld"].notna() & (df["score"] >= df["thrshld"])]
    df = df.rename(columns={{"gene_name": "orf", "KO": "ko", "score": "raw_score"}})
    df = df.merge(ko_to_mnxr, on="ko", how="inner")
    df["source"] = SOURCE
    df["channel"] = "kofam"
    df["intermediate_id"] = df["ko"]
    df["intermediate_name"] = ""
    df["projection_via"] = "kegg.reaction"
    return df[SCHEMA]

# ---- CLEAN lane: dev2 clean_predictions = Query ID \t Predicted EC number \t clean_score
def lane_clean(path, ec_to_mnxr):
    df = pd.read_csv(path, sep="\t")
    df.columns = ["orf", "ec", "clean_score"][:len(df.columns)]
    df["raw_score"] = pd.to_numeric(df["clean_score"], errors="coerce")
    df = df[df["ec"].astype(str).str.match(r"^\d+\.\d+\.\d+\.\d+$", na=False)]
    df = df[df["raw_score"] >= 0.01]                    # CLEAN never abstains: low floor
    df = df.merge(ec_to_mnxr, on="ec", how="inner")
    df["source"] = SOURCE
    df["channel"] = "clean"
    df["intermediate_id"] = df["ec"]
    df["intermediate_name"] = ""
    df["projection_via"] = "ec"
    return df[SCHEMA]

# ---- uniref50 lane: dev2 diamond_uniref50_results = BLAST6 + stitle + bsr (14 col)
def lane_uniref(path, uniprot_to_mnxr):
    df = pd.read_csv(path, sep="\t", header=None, names=fe._BLAST6_BSR_COLS, dtype=str)
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["bsr"] = pd.to_numeric(df["bsr"], errors="coerce")
    df = (df.sort_values(["qseqid", "evalue", "bitscore"], ascending=[True, True, False])
            .drop_duplicates(subset=["qseqid"], keep="first"))
    df["uniprot_accession"] = df["sseqid"].str.replace(r"^UniRef50_", "", regex=True)
    df["intermediate_name"] = df["stitle"].apply(fe._clean_stitle)
    df = df.rename(columns={{"qseqid": "orf", "bsr": "raw_score"}})
    joined = df.merge(uniprot_to_mnxr[["uniprot_accession", "dr_source", "mnxr"]],
                      on="uniprot_accession", how="inner")
    joined["source"] = SOURCE
    joined["channel"] = "uniref50"
    joined["intermediate_id"] = joined["uniprot_accession"]
    joined["projection_via"] = joined["dr_source"]
    return joined[SCHEMA]

# ---- embedding lane: numpy kNN label transfer vs the reference pool ----
# lightweight port of fabfos_embed_transfer.apply_one (cosine top-K distance vote).
def _load_query(parquet, index_csv):
    idx = pd.read_csv(index_csv)
    q = pd.read_parquet(parquet).to_numpy(dtype=np.float32)
    return idx["sequence_id"].to_numpy(), q

def _norm(x):
    n = np.linalg.norm(x, axis=1, keepdims=True)
    return x / np.clip(n, 1e-9, None)

def lane_embed(parquet, index_csv, pool_dir, emb_name, channel, floor):
    pool_idx = pd.read_parquet(os.path.join(pool_dir, "orf_index.parquet"))
    ref = pool_idx[pool_idx["role"] == "reference"].reset_index(drop=True)
    emb = np.load(os.path.join(pool_dir, emb_name), mmap_mode="r")
    ref_emb = _norm(np.asarray(emb[ref["row"].to_numpy()], dtype=np.float32))
    ref_orf = ref["orf"].to_numpy()
    label_lists = [s.split(";") if s else [] for s in ref["mnxr_list"]]
    vocab = sorted({{m for ls in label_lists for m in ls}})
    vidx = {{m: i for i, m in enumerate(vocab)}}
    L = np.zeros((len(ref), len(vocab)), dtype=np.float32)
    for r, ls in enumerate(label_lists):
        for m in ls:
            L[r, vidx[m]] = 1.0
    q_orf, q_emb = _load_query(parquet, index_csv)
    q_emb = _norm(q_emb)
    rows = []
    for s in range(0, len(q_emb), 256):
        sim = q_emb[s:s+256] @ ref_emb.T
        top = np.argpartition(-sim, min(K, sim.shape[1]-1), axis=1)[:, :K]
        for bi in range(sim.shape[0]):
            nn = top[bi]
            vals = np.clip(sim[bi, nn], 0, None)
            tot = vals.sum()
            if tot <= 0:
                continue
            w = vals / tot
            votes = w @ L[nn]
            best = int(np.argmax(sim[bi, nn]))
            for j in np.nonzero(votes >= floor)[0]:
                rows.append((SOURCE, q_orf[s+bi], channel, vocab[j],
                             ref_orf[nn[best]], "", float(votes[j]), "embedding_knn"))
    return pd.DataFrame(rows, columns=SCHEMA)

def main():
    ids = set(orf_ids("{orfs}"))
    ko_to_mnxr = fe.load_ko_to_mnxr("{ko_br}")
    ec_to_mnxr = pd.read_csv("{ec_br}", sep="\t")[["ec", "mnxr"]].drop_duplicates()
    up_to_mnxr = fe.load_uniprot_to_mnxr("{up_br}")
    frames = [
        lane_kofam("{kofam}", ko_to_mnxr),
        lane_clean("{clean}", ec_to_mnxr),
        lane_uniref("{uniref}", up_to_mnxr),
        lane_embed("{pbert_emb}", "{pbert_idx}", "{pool}", "emb_pbert.npy", "pbert", PBERT_FLOOR),
    ]
    gpr = pd.concat(frames, ignore_index=True)
    gpr = gpr[gpr["orf"].isin(ids)]
    gpr.to_parquet("{out}", index=False)
    print(f"[gpr_4lane] wrote {{len(gpr)}} rows ({{gpr['channel'].nunique()}} channels)", flush=True)

main()
'''


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    ikof  = context.Input(kofam)
    icln  = context.Input(clean)
    iuni  = context.Input(uniref)
    ipe   = context.Input(pbert_emb)
    ipi   = context.Input(pbert_idx)
    ikb   = context.Input(ko_br)
    ieb   = context.Input(ec_br)
    iub   = context.Input(up_br)
    ipool = context.Input(pool)
    iev   = context.Input(ev_lib)
    iknn  = context.Input(knn_lib)
    iout  = context.Output(out_gpr)

    driver = DRIVER.format(
        ev_lib=iev.container, knn_lib=iknn.container,
        orfs=iorfs.container, kofam=ikof.container, clean=icln.container,
        uniref=iuni.container, pbert_emb=ipe.container, pbert_idx=ipi.container,
        ko_br=ikb.container, ec_br=ieb.container, up_br=iub.container,
        pool=ipool.container, out=iout.container,
    )
    context.LocalShell("cat > _gpr_4lane.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithContainer(image=image, cmd="python3 _gpr_4lane.py")

    return ExecutionResult(
        manifest=[{out_gpr: iout.local}],
        success=iout.local.exists(),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=orfs,
    resources=Resources(
        cpus=2,
        memory=Size.GB(8),
        duration=Duration(hours=1),
    ),
)
