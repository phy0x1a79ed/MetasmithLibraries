"""GPR mapper (full 7 lanes) -- fold every annotation lane into one gene-attributed
GPR table -> annotation::gpr_table_7lane.

Same shape and contract as gpr_4lane (see that file's header for the 8-col schema
and the "run tool -> mapper" rationale); this is the full-coverage variant that
adds the three lanes the chosen-4 subset omits:

  DeepEC  (EC->MNXR, score-less -> raw_score = NaN)
  EZpred  (enzyme-head level-4 EC->MNXR, softmax score)
  ESM-C   (embedding kNN transfer, like ProteinBERT but the emb_esmc pool)

All seven channels land in one long-format parquet; no cross-lane dedup. Reuses
resources/lib/fabfos_evidence.py (bridge loaders + schema + floors, pandas-only)
and a lightweight numpy port of the fabfos_embed_transfer kNN vote.

GATED: same staged-reference gating as gpr_4lane; planning is type-driven, running
end-to-end is a later milestone.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
orfs      = model.AddRequirement(lib.GetType("sequences::orfs"))
kofam     = model.AddRequirement(lib.GetType("annotation::kofamscan_results"))
clean     = model.AddRequirement(lib.GetType("annotation::clean_predictions"))
deepec    = model.AddRequirement(lib.GetType("annotation::deepec_predictions"))
ezpred    = model.AddRequirement(lib.GetType("annotation::ezpred_predictions"))
uniref    = model.AddRequirement(lib.GetType("annotation::diamond_uniref50_results"))
pbert_emb = model.AddRequirement(lib.GetType("annotation::proteinbert_embeddings"))
pbert_idx = model.AddRequirement(lib.GetType("annotation::proteinbert_index"))
esmc_emb  = model.AddRequirement(lib.GetType("annotation::esm_c_embeddings"))
esmc_idx  = model.AddRequirement(lib.GetType("annotation::esm_c_index"))
bridge    = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
pool      = model.AddRequirement(lib.GetType("ref::reference_label_pool"))
ev_lib    = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
knn_lib   = model.AddRequirement(lib.GetType("lib::fabfos_embed_transfer.py"))
out_gpr   = model.AddProduct(lib.GetType("annotation::gpr_table_7lane"))


DRIVER = r'''
import sys, os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname("{ev_lib}"))
import fabfos_evidence as fe

SCHEMA = fe.SCHEMA_COLS
SOURCE = "fosmid"
K = 30
PBERT_FLOOR = 0.20
ESMC_FLOOR = 0.10
DL_EC_FLOOR = fe.DL_EC_SCORE_FLOOR

def orf_ids(fasta):
    out = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].split()[0])
    return out

def lane_kofam(path, ko_to_mnxr):
    df = pd.read_csv(path)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["thrshld"] = pd.to_numeric(df["thrshld"], errors="coerce")
    df = df[df["score"].notna() & df["thrshld"].notna() & (df["score"] >= df["thrshld"])]
    df = df.rename(columns={{"gene_name": "orf", "KO": "ko", "score": "raw_score"}})
    df = df.merge(ko_to_mnxr, on="ko", how="inner")
    df["source"] = SOURCE; df["channel"] = "kofam"
    df["intermediate_id"] = df["ko"]; df["intermediate_name"] = ""
    df["projection_via"] = "kegg.reaction"
    return df[SCHEMA]

def lane_clean(path, ec_to_mnxr):
    df = pd.read_csv(path, sep="\t")
    df.columns = ["orf", "ec", "clean_score"][:len(df.columns)]
    df["raw_score"] = pd.to_numeric(df["clean_score"], errors="coerce")
    df = df[df["ec"].astype(str).str.match(r"^\d+\.\d+\.\d+\.\d+$", na=False)]
    df = df[df["raw_score"] >= 0.01]
    df = df.merge(ec_to_mnxr, on="ec", how="inner")
    df["source"] = SOURCE; df["channel"] = "clean"
    df["intermediate_id"] = df["ec"]; df["intermediate_name"] = ""
    df["projection_via"] = "ec"
    return df[SCHEMA]

# DeepEC: dev2 deepec_predictions is a TSV; col0 = gene id, col1 = predicted EC.
# score-less -> raw_score = NaN (accepted).
def lane_deepec(path, ec_to_mnxr):
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene, ec = parts[0].strip(), parts[1].strip()
            if not ec or ec.lower().startswith("predicted") or ec.count(".") < 1:
                continue
            rows.append((gene, ec))
    df = pd.DataFrame(rows, columns=["orf", "ec"])
    df = df[df["ec"].str.match(r"^\d+\.\d+\.\d+\.\d+$", na=False)]
    df = df.merge(ec_to_mnxr, on="ec", how="inner")
    df["source"] = SOURCE; df["channel"] = "deepec"
    df["intermediate_id"] = df["ec"]; df["intermediate_name"] = ""
    df["raw_score"] = np.nan; df["projection_via"] = "ec"
    return df[SCHEMA]

# EZpred: dev2 ezpred_predictions = sequence_id, ec_number, score, head_kind.
def lane_ezpred(path, ec_to_mnxr):
    df = pd.read_csv(path)
    df = df[df["head_kind"] == "enzyme"]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df[df["score"] >= DL_EC_FLOOR]
    df = df[df["ec_number"].astype(str).str.match(r"^\d+\.\d+\.\d+\.\d+$", na=False)]
    df = df.merge(ec_to_mnxr, left_on="ec_number", right_on="ec", how="inner")
    df = df.rename(columns={{"sequence_id": "orf", "score": "raw_score"}})
    df["source"] = SOURCE; df["channel"] = "ezpred"
    df["intermediate_id"] = df["ec_number"]; df["intermediate_name"] = ""
    df["projection_via"] = "ec"
    return df[SCHEMA]

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
    joined["source"] = SOURCE; joined["channel"] = "uniref50"
    joined["intermediate_id"] = joined["uniprot_accession"]
    joined["projection_via"] = joined["dr_source"]
    return joined[SCHEMA]

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
    idx = pd.read_csv(index_csv)
    q_orf = idx["sequence_id"].to_numpy()
    q_emb = _norm(pd.read_parquet(parquet).to_numpy(dtype=np.float32))
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

def load_bridge(path):
    """One table, three id spaces. Sliced by id_source into the per-lane frames the
    lane functions expect. The spaces share no ids, so the slice is exact."""
    b = pd.read_parquet(path, columns=["id", "id_source", "mnxr"])
    def slice_as(src, name):
        s = b[b["id_source"] == src][["id", "mnxr"]].drop_duplicates()
        return s.rename(columns={{"id": name}})
    return (slice_as("ko", "ko"),
            slice_as("ec", "ec"),
            slice_as("uniprot", "uniprot_accession").assign(dr_source="rhea"))

def main():
    ids = set(orf_ids("{orfs}"))
    ko_to_mnxr, ec_to_mnxr, up_to_mnxr = load_bridge("{bridge}")
    frames = [
        lane_kofam("{kofam}", ko_to_mnxr),
        lane_clean("{clean}", ec_to_mnxr),
        lane_deepec("{deepec}", ec_to_mnxr),
        lane_ezpred("{ezpred}", ec_to_mnxr),
        lane_uniref("{uniref}", up_to_mnxr),
        lane_embed("{pbert_emb}", "{pbert_idx}", "{pool}", "emb_pbert.npy", "pbert", PBERT_FLOOR),
        lane_embed("{esmc_emb}", "{esmc_idx}", "{pool}", "emb_esmc.npy", "esmc", ESMC_FLOOR),
    ]
    gpr = pd.concat(frames, ignore_index=True)
    gpr = gpr[gpr["orf"].isin(ids)]
    gpr.to_parquet("{out}", index=False)
    print(f"[gpr_7lane] wrote {{len(gpr)}} rows ({{gpr['channel'].nunique()}} channels)", flush=True)

main()
'''


def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    ikof  = context.Input(kofam)
    icln  = context.Input(clean)
    idec  = context.Input(deepec)
    iez   = context.Input(ezpred)
    iuni  = context.Input(uniref)
    ipe   = context.Input(pbert_emb)
    ipi   = context.Input(pbert_idx)
    iee   = context.Input(esmc_emb)
    iei   = context.Input(esmc_idx)
    ibr   = context.Input(bridge)
    ipool = context.Input(pool)
    iev   = context.Input(ev_lib)
    iknn  = context.Input(knn_lib)
    iout  = context.Output(out_gpr)

    driver = DRIVER.format(
        ev_lib=iev.container, knn_lib=iknn.container,
        orfs=iorfs.container, kofam=ikof.container, clean=icln.container,
        deepec=idec.container, ezpred=iez.container, uniref=iuni.container,
        pbert_emb=ipe.container, pbert_idx=ipi.container,
        esmc_emb=iee.container, esmc_idx=iei.container,
        bridge=ibr.container,
        pool=ipool.container, out=iout.container,
    )
    context.LocalShell("cat > _gpr_7lane.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithContainer(image=image, cmd="python3 _gpr_7lane.py")

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
