"""ECSPr Layer-2 evidence: MNXR bridges + 4-lane readers + evidence-table compiler.

Faithful port of the scadc `03_layer2_evidence/` pipeline into a single
path-parameterized CLI so a metasmith transform can call it (no hardcoded scadc
paths). Ported verbatim (method-preserving) from:
  - _lanes.py                     (bridge loaders + the 4 lane readers, 8-col schema)
  - 00_build_uniprot_to_mnxr.py   (UniProt->MNXR Rhea DR bridge)
  - 11_build_evidence_table_dlec.py (the dl_ec 4-channel compiler)

Every lane reader emits the unified 8-column schema:
    source, orf, channel, mnxr, intermediate_id, intermediate_name, raw_score, projection_via

The compiler folds the four fresh-annotation lanes into one evidence table (the
canonical `dl_ec` variant); ECSPr consumes it as functional_annotation::evidence_table.

SUBCOMMANDS
  build-uniprot-bridge  reac_xref + rhea2uniprot{,_trembl}  -> uniprot_to_mnxr.parquet
  build-ec-bridge       reac_prop.tsv                       -> ec_to_mnxr.tsv
  compile               per-lane annotator outputs + bridges -> evidence_table.parquet

ko_to_mnxr is a reused reference table (staged input), not built here.

GATES: the lane INPUTS are produced fresh upstream by the heavy annotator
transforms (kofamscan, EZpred/ESM-C dl_ec, DIAMOND-vs-UniRef50, embed-transfer).
This module only needs pandas/pyarrow + the MetaNetX/Rhea reference tables (present
on disk); it recomputes NOTHING that must be reused (the metaG null lives elsewhere).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

SCHEMA_COLS = [
    "source", "orf", "channel", "mnxr",
    "intermediate_id", "intermediate_name", "raw_score", "projection_via",
]

# EZpred `dl_ec` lane keeps level-4 ECs whose softmax score clears this floor.
DL_EC_SCORE_FLOOR = 0.3
# embed-transfer lane keeps pbert_transfer calls whose kNN vote fraction clears this.
EMBED_SCORE_FLOOR = 0.2


# =====================================================================
# bridge loaders
# =====================================================================

def load_ko_to_mnxr(path: Path) -> pd.DataFrame:
    """KO -> MNXR fan-out table (reused reference cache). Multi-row per KO is normal."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"ko": "ko", "mnx_r": "mnxr"})
    df = df[["ko", "mnxr"]].drop_duplicates()
    return df


def load_ec_to_mnxr(path: Path) -> pd.DataFrame:
    """Parse reac_prop.tsv classifs (col 4) -> long-format ec -> mnxr fan-out.

    classifs is `;`-separated; each token is a raw EC number (no `ec:` prefix).
    Partial ECs like "1.2.1" are kept verbatim.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            mnxr, classifs = parts[0], parts[3]
            if not mnxr.startswith("MNXR") or not classifs:
                continue
            for ec in classifs.split(";"):
                ec = ec.strip()
                if not ec:
                    continue
                rows.append((ec, mnxr))
    return pd.DataFrame(rows, columns=["ec", "mnxr"]).drop_duplicates()


def load_uniprot_to_mnxr(path: Path) -> pd.DataFrame:
    """The cache produced by build-uniprot-bridge."""
    return pd.read_parquet(path)


# =====================================================================
# helpers
# =====================================================================

_BOILERPLATE_RE = re.compile(r"\s+(n=\d+|Tax=.+?|RepID=\S+)(?=\s|$)")


def _clean_stitle(stitle: str) -> str:
    """Strip UniRef50 boilerplate from a stitle to get a readable name."""
    if not isinstance(stitle, str):
        return ""
    s = stitle
    if s.startswith("UniRef50_"):
        head, _, rest = s.partition(" ")
        s = rest
    return _BOILERPLATE_RE.sub("", s).strip()


# =====================================================================
# lane readers (unified 8-col schema)
# =====================================================================

def read_kofam(path, source: str, ko_to_mnxr: pd.DataFrame) -> pd.DataFrame:
    """Load *.kofam.csv, keep above-threshold hits, project KO -> MNXR."""
    if path is None or not Path(path).exists():
        return pd.DataFrame(columns=SCHEMA_COLS)
    df = pd.read_csv(path)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["hmm_threshold"] = pd.to_numeric(df["hmm_threshold"], errors="coerce")
    df = df[df["score"].notna() & df["hmm_threshold"].notna()]
    df = df[df["score"] >= df["hmm_threshold"]]
    if "fosmid" in df.columns:
        df["orf_id"] = df["fosmid"].astype(str) + "_" + df["orf"].astype(str)
    elif "contig" in df.columns:
        df["orf_id"] = df["contig"].astype(str) + "_" + df["orf"].astype(str)
    else:
        df["orf_id"] = df["orf"].astype(str)
    df = df[["orf_id", "ko", "score", "description"]].rename(
        columns={"orf_id": "orf", "score": "raw_score", "description": "intermediate_name"}
    )
    df = df.merge(ko_to_mnxr, on="ko", how="inner")
    df["source"] = source
    df["channel"] = "kofam"
    df["intermediate_id"] = df["ko"]
    df["projection_via"] = "kegg.reaction"
    return df[SCHEMA_COLS]


def read_dl_ec(path, source: str, ec_to_mnxr: pd.DataFrame) -> pd.DataFrame:
    """Load EZpred (ESM-C 600M DL-only) EC predictions, project EC -> MNXR.

    parquet columns: sequence_id, ec_number, score, head_kind. Keeps enzyme-head,
    level-4 ECs (x.x.x.x) clearing DL_EC_SCORE_FLOOR; raw_score carries the
    EZpred confidence. Filters pushed to pyarrow so huge parquets never fully load.
    """
    if path is None or not Path(path).exists():
        return pd.DataFrame(columns=SCHEMA_COLS)
    df = pd.read_parquet(
        path,
        columns=["sequence_id", "ec_number", "score"],
        filters=[("head_kind", "==", "enzyme"), ("score", ">=", DL_EC_SCORE_FLOOR)],
    )
    df = df[df["ec_number"].astype(str).str.match(r"^\d+\.\d+\.\d+\.\d+$", na=False)]
    df = df.merge(ec_to_mnxr, left_on="ec_number", right_on="ec", how="inner")
    df["source"] = source
    df["channel"] = "dl_ec"
    df["intermediate_id"] = df["ec_number"]
    df["intermediate_name"] = ""
    df = df.rename(columns={"sequence_id": "orf", "score": "raw_score"})
    df["orf"] = df["orf"].str.replace(r"-(\d+)$", r"_\1", regex=True)
    df["projection_via"] = "ec"
    return df[SCHEMA_COLS]


_BLAST6_BSR_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle", "bsr",
]


def read_uniref50(path, source: str, uniprot_to_mnxr: pd.DataFrame) -> pd.DataFrame:
    """Load DIAMOND BLAST6+stitle+BSR; best-hit per ORF; project UniProt -> MNXR."""
    if path is None or not Path(path).exists():
        return pd.DataFrame(columns=SCHEMA_COLS)
    df = pd.read_csv(path, sep="\t", header=None, names=_BLAST6_BSR_COLS, dtype=str)
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["bsr"] = pd.to_numeric(df["bsr"], errors="coerce")
    df = (df.sort_values(["qseqid", "evalue", "bitscore"], ascending=[True, True, False])
            .drop_duplicates(subset=["qseqid"], keep="first"))
    df["uniprot_accession"] = df["sseqid"].str.replace(r"^UniRef50_", "", regex=True)
    df["intermediate_name"] = df["stitle"].apply(_clean_stitle)
    df = df[["qseqid", "uniprot_accession", "intermediate_name", "bsr"]].rename(
        columns={"qseqid": "orf", "bsr": "raw_score"})
    df["orf"] = df["orf"].str.replace(r"-(\d+)$", r"_\1", regex=True)
    joined = df.merge(uniprot_to_mnxr[["uniprot_accession", "dr_source", "mnxr"]],
                      on="uniprot_accession", how="inner")
    joined["source"] = source
    joined["channel"] = "uniref50_dr"
    joined["intermediate_id"] = joined["uniprot_accession"]
    joined["projection_via"] = joined["dr_source"]
    return joined[SCHEMA_COLS]


def read_embed_transfer(path, source: str, _bridge=None) -> pd.DataFrame:
    """Load the embedding-transfer (lane 4) candidate table, ProteinBERT channel.

    The candidate parquet already emits the 8-col schema (+ audit cols); filter to
    the pbert_transfer channel above the vote floor. Label transfer *is* the
    projection (mnxr already present, projection_via='embedding_knn').
    """
    if path is None or not Path(path).exists():
        return pd.DataFrame(columns=SCHEMA_COLS)
    df = pd.read_parquet(path)
    df = df[(df["channel"] == "pbert_transfer") & (df["raw_score"] >= EMBED_SCORE_FLOOR)].copy()
    df["source"] = source
    return df[SCHEMA_COLS]


# =====================================================================
# build-uniprot-bridge  (port of 00_build_uniprot_to_mnxr.py)
# =====================================================================

def _load_rhea_to_mnxr(reac_xref: Path) -> pd.DataFrame:
    rows = []
    with open(reac_xref) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("EMPTY"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, mnxr = parts[0], parts[1]
            if not src.startswith("rhea:") or not mnxr.startswith("MNXR"):
                continue
            rows.append((int(src[len("rhea:"):]), mnxr))
    return pd.DataFrame(rows, columns=["rhea_id", "mnxr"])


def _load_rhea2uniprot(path: Path, evidence_quality: str) -> pd.DataFrame:
    df = pd.read_csv(
        path, sep="\t",
        dtype={"RHEA_ID": "int64", "DIRECTION": "string", "MASTER_ID": "int64", "ID": "string"},
    )
    df = df.rename(columns={"RHEA_ID": "rhea_id", "DIRECTION": "direction", "ID": "uniprot_accession"})
    df = df[["rhea_id", "direction", "uniprot_accession"]]
    df["evidence_quality"] = evidence_quality
    return df


def build_uniprot_bridge(reac_xref: Path, rhea_swiss: Path, rhea_trembl: Path, out: Path):
    print("[bridge] rhea -> MNXR from reac_xref...", flush=True)
    rhea_mnxr = _load_rhea_to_mnxr(reac_xref)
    print(f"         {len(rhea_mnxr):,} rows ({rhea_mnxr['mnxr'].nunique():,} MNXRs)", flush=True)
    frames = [_load_rhea2uniprot(rhea_swiss, "reviewed")]
    if rhea_trembl is not None and Path(rhea_trembl).exists():
        print("[bridge] rhea2uniprot TrEMBL (big)...", flush=True)
        frames.append(_load_rhea2uniprot(rhea_trembl, "unreviewed"))
    all_u = pd.concat(frames, ignore_index=True)
    joined = all_u.merge(rhea_mnxr, on="rhea_id", how="inner")
    out_df = pd.DataFrame({
        "uniprot_accession": joined["uniprot_accession"],
        "protein_name": "", "gene_name": "",
        "dr_source": "rhea",
        "external_id": "rhea:" + joined["rhea_id"].astype(str),
        "mnxr": joined["mnxr"],
        "evidence_quality": joined["evidence_quality"],
        "direction": joined["direction"],
    })
    out_df = (out_df.sort_values(["uniprot_accession", "mnxr", "evidence_quality"])
                    .drop_duplicates(subset=["uniprot_accession", "external_id", "mnxr"], keep="first"))
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out, index=False)
    print(f"[bridge] wrote {len(out_df):,} rows ({out_df['uniprot_accession'].nunique():,} UniProts) -> {out}", flush=True)


# =====================================================================
# compile  (port of 11_build_evidence_table_dlec.py, single-source capable)
# =====================================================================

def compile_evidence(source, kofam, dl_ec, uniref50, embed,
                     ko_to_mnxr_path, ec_to_mnxr_path, uniprot_to_mnxr_path, out):
    print("[compile] loading bridges...", flush=True)
    ko_to_mnxr = load_ko_to_mnxr(ko_to_mnxr_path) if ko_to_mnxr_path else pd.DataFrame(columns=["ko", "mnxr"])
    # ec_to_mnxr here is the PRE-BUILT 2-col bridge (from build-ec-bridge), not
    # reac_prop -- read it plainly rather than re-running the reac_prop parser.
    if ec_to_mnxr_path:
        ec_to_mnxr = pd.read_csv(ec_to_mnxr_path, sep="\t")[["ec", "mnxr"]].drop_duplicates()
    else:
        ec_to_mnxr = pd.DataFrame(columns=["ec", "mnxr"])
    if uniprot_to_mnxr_path and Path(uniprot_to_mnxr_path).exists():
        uniprot_to_mnxr = load_uniprot_to_mnxr(uniprot_to_mnxr_path)
    else:
        print("[compile] uniprot_to_mnxr absent -> uniref50 lane skipped", flush=True)
        uniprot_to_mnxr = pd.DataFrame(columns=["uniprot_accession", "dr_source", "mnxr"])

    lanes = [
        ("kofam", read_kofam, kofam, ko_to_mnxr),
        ("dl_ec", read_dl_ec, dl_ec, ec_to_mnxr),
        ("uniref50", read_uniref50, uniref50, uniprot_to_mnxr),
        ("embed", read_embed_transfer, embed, None),
    ]
    frames = []
    for name, reader, path, bridge in lanes:
        if not path:
            continue
        f = reader(path, source, bridge)
        print(f"[compile] {name}: {len(f):,} rows ({f['orf'].nunique():,} ORFs, {f['mnxr'].nunique():,} MNXRs)", flush=True)
        frames.append(f)
    if not frames:
        raise SystemExit("no lane inputs provided; nothing to compile")
    out_df = pd.concat(frames, ignore_index=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(out, index=False)
    print(f"[compile] wrote {len(out_df):,} rows -> {out}", flush=True)


# =====================================================================
# CLI
# =====================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build-uniprot-bridge", help="UniProt->MNXR via Rhea DR")
    b.add_argument("--reac-xref", type=Path, required=True)
    b.add_argument("--rhea-swiss", type=Path, required=True)
    b.add_argument("--rhea-trembl", type=Path, default=None)
    b.add_argument("--out", type=Path, required=True)

    e = sub.add_parser("build-ec-bridge", help="EC(level-4)->MNXR from reac_prop classifs")
    e.add_argument("--reac-prop", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)

    c = sub.add_parser("compile", help="fold lanes -> evidence_table.parquet")
    c.add_argument("--source", default="fosmid")
    c.add_argument("--kofam", type=Path, default=None)
    c.add_argument("--dl-ec", type=Path, default=None)
    c.add_argument("--uniref50", type=Path, default=None)
    c.add_argument("--embed", type=Path, default=None)
    c.add_argument("--ko-to-mnxr", type=Path, default=None)
    c.add_argument("--ec-to-mnxr", type=Path, default=None)
    c.add_argument("--uniprot-to-mnxr", type=Path, default=None)
    c.add_argument("--out", type=Path, required=True)

    a = ap.parse_args()
    if a.cmd == "build-uniprot-bridge":
        build_uniprot_bridge(a.reac_xref, a.rhea_swiss, a.rhea_trembl, a.out)
    elif a.cmd == "build-ec-bridge":
        df = load_ec_to_mnxr(a.reac_prop)
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(a.out, sep="\t", index=False)
        print(f"[bridge] wrote {len(df):,} ec->mnxr rows ({df['ec'].nunique():,} ECs) -> {a.out}", flush=True)
    elif a.cmd == "compile":
        compile_evidence(a.source, a.kofam, a.dl_ec, a.uniref50, a.embed,
                         a.ko_to_mnxr, a.ec_to_mnxr, a.uniprot_to_mnxr, a.out)


if __name__ == "__main__":
    main()
