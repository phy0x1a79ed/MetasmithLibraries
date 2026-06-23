from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
img_blast = model.AddRequirement(lib.GetType("containers::blast.oci"))
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
asm       = model.AddRequirement(lib.GetType("sequences::assembly"))
table     = model.AddRequirement(lib.GetType("fosmids::end_sequences_table"))
fwd_ends  = model.AddRequirement(lib.GetType("fosmids::end_sequences_forward"), parents={table})
rev_ends  = model.AddRequirement(lib.GetType("fosmids::end_sequences_reverse"), parents={table})
out_fa    = model.AddProduct(lib.GetType("fosmids::scaffolds"))
out_rep   = model.AddProduct(lib.GetType("fosmids::scaffold_report"))

def protocol(context: ExecutionContext):
    iasm   = context.Input(asm)
    itable = context.Input(table)
    ifwd   = context.Input(fwd_ends)
    irev   = context.Input(rev_ends)
    ofa    = context.Output(out_fa)
    orep   = context.Output(out_rep)

    exp_length       = context.params.get("exp_length", 40000)
    exp_length_range = context.params.get("exp_length_range", 5000)
    min_contig_len   = context.params.get("min_contig_length", 500)
    gap_str          = context.params.get("gap_str", "N" * 100)
    pident_thresh    = context.params.get("pident", 90)
    threads          = context.params.get("cpus")
    threads_arg      = "" if threads is None else f"-num_threads {threads}"

    # ---------------------------------------------------------------
    # Stage 1: prepare contigs + ends in workspace
    prep_script = "prep.py"
    with open(prep_script, "w") as f:
        f.write(f"""\
import sys, json, math
from pathlib import Path
from Bio import SeqIO

asm_path = sys.argv[1]
fwd_path = sys.argv[2]
rev_path = sys.argv[3]
table_path = sys.argv[4]

with open(table_path) as j:
    table = json.load(j)
insert_ids = list(table.get("insert_ids", []))
ends_facing = bool(table.get("ends_facing", False))

# normalize reverse ends so both fwd and rev face the insert interior
prepared_rev = "prepared_rev.fna"
if ends_facing:
    with open(prepared_rev, "w") as out:
        for rec in SeqIO.parse(rev_path, "fasta"):
            out.write(f">{{rec.id}}\\n{{str(rec.seq)}}\\n")
else:
    with open(prepared_rev, "w") as out:
        for rec in SeqIO.parse(rev_path, "fasta"):
            out.write(f">{{rec.id}}\\n{{str(rec.seq.reverse_complement())}}\\n")

# copy fwd as-is (with normalized header)
prepared_fwd = "prepared_fwd.fna"
with open(prepared_fwd, "w") as out:
    for rec in SeqIO.parse(fwd_path, "fasta"):
        out.write(f">{{rec.id}}\\n{{str(rec.seq)}}\\n")

# build all_contigs with stable rekeyed IDs (C0001 ...)
contig_records = [r for r in SeqIO.parse(asm_path, "fasta") if len(r.seq) >= {min_contig_len}]
n = len(contig_records)
assert n > 0, f"no contigs >= {min_contig_len}bp"
digits = max(1, int(math.ceil(math.log10(n + 1))))
id_map = {{}}
contig_lengths = {{}}
contig_seqs = {{}}
with open("all_contigs.fa", "w") as out:
    for i, rec in enumerate(contig_records, start=1):
        k = f"C{{i:0{{digits}}d}}"
        id_map[k] = rec.id
        contig_lengths[k] = len(rec.seq)
        contig_seqs[k] = str(rec.seq)
        out.write(f">{{k}} length={{len(rec.seq)}} | {{rec.description}}\\n{{str(rec.seq)}}\\n")

with open("prep_state.json", "w") as j:
    json.dump(dict(
        insert_ids=insert_ids,
        id_map=id_map,
        contig_lengths=contig_lengths,
    ), j)

# emit contig sequences as a separate fasta-like JSON for fast lookup downstream
with open("contig_seqs.json", "w") as j:
    json.dump(contig_seqs, j)

print(f"contigs={{n}} fwd_ends_path={{prepared_fwd}} rev_ends_path={{prepared_rev}}")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {prep_script} {iasm.container} {ifwd.container} {irev.container} {itable.container}",
    )

    # ---------------------------------------------------------------
    # Stage 2: blast fwd & rev ends against contigs
    context.ExecWithContainer(
        image=img_blast,
        cmd=f"""\
            mkdir -p blast_db
            makeblastdb -dbtype nucl -in all_contigs.fa -out blast_db/all_contigs >makeblastdb.log 2>&1
            blastn -evalue 999 -perc_identity {pident_thresh} {threads_arg} \
                -query prepared_fwd.fna \
                -db blast_db/all_contigs \
                -outfmt "6 qseqid sseqid nident length qlen slen qstart qend sstart send" \
                -out hits_forward.tsv >blast_fwd.log 2>&1
            blastn -evalue 999 -perc_identity {pident_thresh} {threads_arg} \
                -query prepared_rev.fna \
                -db blast_db/all_contigs \
                -outfmt "6 qseqid sseqid nident length qlen slen qstart qend sstart send" \
                -out hits_reverse.tsv >blast_rev.log 2>&1
        """,
    )

    # ---------------------------------------------------------------
    # Stage 3: parse hits, build candidate seqs, write scaffolding queries
    parse_script = "parse_hits.py"
    with open(parse_script, "w") as f:
        f.write(f"""\
import json
from enum import Enum
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq as RawSeq

PIDENT_THRESH = {pident_thresh}
EXP_LEN = {exp_length}
EXP_RANGE = {exp_length_range}
GAP_STR = {gap_str!r}

with open("prep_state.json") as j:
    state = json.load(j)
insert_ids = state["insert_ids"]
id_map = state["id_map"]
contig_lengths = state["contig_lengths"]
with open("contig_seqs.json") as j:
    contig_seqs_raw = json.load(j)

# end seqs (use Seq objects so reverse_complement works in the next phase if needed)
def _load_fa(path):
    out = {{}}
    for rec in SeqIO.parse(path, "fasta"):
        out[rec.id] = str(rec.seq)
    return out
fwd_ends = _load_fa("prepared_fwd.fna")
rev_ends = _load_fa("prepared_rev.fna")

class QUALITY(float, Enum):
    FULL_MATCH = 10
    FULL_SCAFFOLD = 7
    GAPPED_SCAFFOLD = 5
    FORWARD_END_ONLY = 3.001
    REVERSE_END_ONLY = 3
    NO_HITS = 0

COLS = "qseqid sseqid nident length qlen slen qstart qend sstart send".split()

def _load_hits(path):
    try:
        df = pd.read_csv(path, sep="\\t", header=None, names=COLS)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=COLS)
    return df

class Hit:
    def __init__(self, row, qseqs, sseqs):
        self.qseqid = row.qseqid
        self.sseqid = row.sseqid
        self.length = int(row.length)
        self.percent_query_mapped = (row.nident / row.qlen) * 100
        qstart, qend = int(row.qstart), int(row.qend)
        sstart, send = int(row.sstart), int(row.send)
        slen = int(row.slen)
        # zero-index, open
        qstart -= 1
        if sstart < send:
            sstart -= 1
        else:
            send -= 1
        if qstart > qend or sstart > send:
            s, e = sstart + qstart, 0
        else:
            s, e = sstart - qstart, slen
        self._subject_trunc_loc = s, e
        self.query_seq = qseqs[self.qseqid]
        self.subject_seq = sseqs[self.sseqid]
        self.subject_len = slen

    @staticmethod
    def _cut(seq_str, start, end):
        if start > end:
            rc = str(RawSeq(seq_str).reverse_complement())
            return rc[len(rc) - start: len(rc) - end]
        return seq_str[start:end]

    def _edges(self, s, e):
        L = self.subject_len
        if s < e:
            lo = abs(min(s, 0)); ro = max(L, e) - L
            return s + lo, e - ro, lo, ro
        lo = max(L, s) - L; ro = abs(min(e, 0))
        return s - lo, e + ro, lo, ro

    def AlignSubjectToQuery(self):
        s, e, lo, ro = self._edges(*self._subject_trunc_loc)
        left  = self._cut(self.query_seq, 0, lo) if lo > 0 else ""
        right = self._cut(self.query_seq, len(self.query_seq) - ro, len(self.query_seq)) if ro > 0 else ""
        mid = self._cut(self.subject_seq, s, e)
        return left + mid + right

    def SeqFromPairedHit(self, other):
        fs, fe = self._subject_trunc_loc
        rs, re = other._subject_trunc_loc
        if (fs > fe) == (rs > re):
            return ""
        if fs < fe and fs > rs: return ""
        if rs < re and rs > fs: return ""
        s, e, lo, ro = self._edges(fs, rs)
        left  = self._cut(self.query_seq, 0, lo) if lo > 0 else ""
        right_rc = self._cut(other.query_seq, ro, 0) if ro > 0 else ""
        return left + self._cut(self.subject_seq, s, e) + right_rc

# bucket hits by insert_id then by contig
def _bucket(df, end_dict):
    out = {{}}
    for _, row in df.iterrows():
        h = Hit(row, end_dict, contig_seqs_raw)
        if h.percent_query_mapped < PIDENT_THRESH:
            continue
        out.setdefault(h.qseqid, []).append(h)
    return out

fwd_hits = _bucket(_load_hits("hits_forward.tsv"), fwd_ends)
rev_hits = _bucket(_load_hits("hits_reverse.tsv"), rev_ends)

mapped_candidates = {{}}  # insert_id -> list[(quality, contig_id, seq)]
to_scaffold = []

def _longest(lst):
    if not lst: return None
    return max(lst, key=lambda t: len(t[2]))

for insert_id in insert_ids:
    fhits = {{}}
    for h in fwd_hits.get(insert_id, []):
        fhits.setdefault(h.sseqid, []).append(h)
    rhits = {{}}
    for h in rev_hits.get(insert_id, []):
        rhits.setdefault(h.sseqid, []).append(h)

    candidates = []
    paired = fhits.keys() & rhits.keys()
    if paired:
        gen = []
        for cid in paired:
            for fh in fhits[cid]:
                for rh in rhits[cid]:
                    seq = fh.SeqFromPairedHit(rh)
                    if seq:
                        gen.append((QUALITY.FULL_MATCH.value, cid, seq))
        best = _longest(gen)
        if best is not None:
            candidates.append(best)

    ufhits = {{k: v for k, v in fhits.items() if k not in rhits}}
    urhits = {{k: v for k, v in rhits.items() if k not in fhits}}
    if ufhits and urhits:
        to_scaffold.append((insert_id, list(ufhits.items()), list(urhits.items())))

    if fhits:
        gen = [(QUALITY.FORWARD_END_ONLY.value, cid, h.AlignSubjectToQuery())
               for cid, hs in fhits.items() for h in hs]
        best = _longest(gen)
        if best is not None:
            candidates.append(best)

    if rhits:
        gen = [(QUALITY.REVERSE_END_ONLY.value, cid, h.AlignSubjectToQuery())
               for cid, hs in rhits.items() for h in hs]
        best = _longest(gen)
        if best is not None:
            candidates.append(best)

    mapped_candidates[insert_id] = candidates

# write fwd/rev "oriented" scaffolding fastas
fwd_sf_seqs = {{}}
rev_sf_seqs = {{}}
SEP = "__"
with open("scaffolding_fwds.fa", "w") as fwf, open("scaffolding_revs.fa", "w") as rvf:
    for insert_id, fl, rl in to_scaffold:
        for cid, hs in fl:
            for h in hs:
                k = f"{{insert_id}}{{SEP}}{{cid}}"
                seq = h.AlignSubjectToQuery()
                fwd_sf_seqs[k] = seq
                fwf.write(f">{{k}}\\n{{seq}}\\n")
        for cid, hs in rl:
            for h in hs:
                k = f"{{insert_id}}{{SEP}}{{cid}}"
                seq = h.AlignSubjectToQuery()
                rev_sf_seqs[k] = seq
                rvf.write(f">{{k}}\\n{{seq}}\\n")

with open("phase2_state.json", "w") as j:
    json.dump(dict(
        mapped=[
            (iid, [[q, cid, seq] for q, cid, seq in cands])
            for iid, cands in mapped_candidates.items()
        ],
        to_scaffold=[
            (iid, [(cid, [h.qseqid for h in hs]) for cid, hs in fl],
                  [(cid, [h.qseqid for h in hs]) for cid, hs in rl])
            for iid, fl, rl in to_scaffold
        ],
    ), j)

print(f"to_scaffold_count={{len(to_scaffold)}}")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {parse_script}",
    )

    # ---------------------------------------------------------------
    # Stage 4: blast fwd-scaffolding-fastas against rev-scaffolding-fastas
    context.ExecWithContainer(
        image=img_blast,
        cmd=f"""\
            if [ -s scaffolding_fwds.fa ] && [ -s scaffolding_revs.fa ]; then
                mkdir -p blast_db_scaffold
                makeblastdb -dbtype nucl -in scaffolding_revs.fa -out blast_db_scaffold/revs >makeblastdb_sf.log 2>&1
                blastn -perc_identity {pident_thresh} {threads_arg} \
                    -query scaffolding_fwds.fa \
                    -db blast_db_scaffold/revs \
                    -outfmt "6 qseqid sseqid nident length qlen slen qstart qend sstart send" \
                    -out scaffolding_blast.tsv >blast_sf.log 2>&1
            else
                : >scaffolding_blast.tsv
            fi
        """,
    )

    # ---------------------------------------------------------------
    # Stage 5: resolve & write outputs
    resolve_script = "resolve.py"
    with open(resolve_script, "w") as f:
        f.write(f"""\
import json, math
import pandas as pd
import numpy as np
from enum import Enum
from Bio.Seq import Seq as RawSeq

EXP_LEN = {exp_length}
EXP_RANGE = {exp_length_range}
GAP_STR = {gap_str!r}
PIDENT_THRESH = {pident_thresh}
MIN_SCAFFOLD_OVERLAP = 150

class QUALITY(float, Enum):
    FULL_MATCH = 10
    FULL_SCAFFOLD = 7
    GAPPED_SCAFFOLD = 5
    FORWARD_END_ONLY = 3.001
    REVERSE_END_ONLY = 3
    NO_HITS = 0

with open("prep_state.json") as j:
    state = json.load(j)
id_map = state["id_map"]
insert_ids = state["insert_ids"]

with open("phase2_state.json") as j:
    p2 = json.load(j)
mapped = {{iid: [tuple(c) for c in cands] for iid, cands in p2["mapped"]}}
to_scaffold = p2["to_scaffold"]

# load scaffolding hits
COLS = "qseqid sseqid nident length qlen slen qstart qend sstart send".split()
try:
    sf_df = pd.read_csv("scaffolding_blast.tsv", sep="\\t", header=None, names=COLS)
except (pd.errors.EmptyDataError, FileNotFoundError):
    sf_df = pd.DataFrame(columns=COLS)

# load oriented scaffolding seqs
def _load_fa(path):
    d = {{}}
    cur = None
    try:
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    d[cur] = []
                else:
                    d[cur].append(line.strip())
    except FileNotFoundError:
        return {{}}
    return {{k: "".join(v) for k, v in d.items()}}

fwd_sf = _load_fa("scaffolding_fwds.fa")
rev_sf = _load_fa("scaffolding_revs.fa")

SEP = "__"

def _split(k):
    p = k.split(SEP)
    if len(p) != 2: return None
    return p[0], p[1]

# union for full scaffolds with overlap
def _union(row):
    qseq = fwd_sf.get(row.qseqid, "")
    sseq = rev_sf.get(row.sseqid, "")
    if not qseq or not sseq: return ""
    qstart, qend = int(row.qstart) - 1, int(row.qend)
    sstart, send = int(row.sstart), int(row.send)
    if sstart < send:
        sstart -= 1
    else:
        send -= 1
    if not (qstart > qend or sstart > send):
        return ""
    qa, qb = 0, qend
    sa, sb = sstart, 0
    overlap = int(row.length)
    qcut = qseq[qa:qb] if qa < qb else str(RawSeq(qseq).reverse_complement())[len(qseq) - qa: len(qseq) - qb]
    if sa > sb:
        rc = str(RawSeq(sseq).reverse_complement())
        scut = rc[len(rc) - sa: len(rc) - sb]
    else:
        scut = sseq[sa:sb]
    return qcut + scut[overlap:]

scaffolding_hits = {{}}
for _, row in sf_df.iterrows():
    fwd_id = _split(row.qseqid)
    rev_id = _split(row.sseqid)
    if fwd_id is None or rev_id is None: continue
    if fwd_id[0] != rev_id[0]: continue
    if fwd_id[1] == rev_id[1]: continue
    if int(row.length) < MIN_SCAFFOLD_OVERLAP: continue
    scaffolding_hits.setdefault(fwd_id[0], []).append((row, fwd_id[1], rev_id[1]))

# build full-scaffold and gapped-scaffold candidates
def _longest(lst):
    if not lst: return None
    return max(lst, key=lambda t: len(t[2]))

for insert_id, fl, rl in to_scaffold:
    new = []
    if insert_id in scaffolding_hits:
        gen = [(QUALITY.FULL_SCAFFOLD.value, f"{{f}};{{r}}", _union(row))
               for row, f, r in scaffolding_hits[insert_id]]
        gen = [t for t in gen if t[2]]
        best = _longest(gen)
        if best is not None:
            new.append(best)
    # gapped scaffolds: every fwd hit (oriented) + GAP + rev hit (rc)
    for cid_f, _ in fl:
        fseq = fwd_sf.get(f"{{insert_id}}{{SEP}}{{cid_f}}", "")
        if not fseq: continue
        for cid_r, _ in rl:
            rseq = rev_sf.get(f"{{insert_id}}{{SEP}}{{cid_r}}", "")
            if not rseq: continue
            joined = fseq + GAP_STR + str(RawSeq(rseq).reverse_complement())
            new.append((QUALITY.GAPPED_SCAFFOLD.value, f"{{cid_f}};{{cid_r}}", joined))
    mapped[insert_id] = mapped.get(insert_id, []) + new

# score & resolve
def _gauss(x):
    if EXP_RANGE == 0: return 0
    k = 1.75
    n = (x - EXP_LEN) / (k * EXP_RANGE)
    return float(np.exp(-(n ** 2)))

def _score(q_val, length):
    qfrac = _gauss(length)
    qscore = q_val / QUALITY.FULL_MATCH.value
    deviation = abs(EXP_LEN - length)
    lscore = max(0.0, 1 - (deviation / EXP_LEN))
    return (qscore * qfrac) + (lscore * (1 - qfrac))

def _name(val):
    for q in QUALITY:
        if q.value == val: return q.name.lower()
    return "unknown"

import sys
mapped_path = sys.argv[1]
report_path = sys.argv[2]

rows = []
with open(mapped_path, "w") as fa:
    for insert_id in sorted(insert_ids):
        cands = mapped.get(insert_id, [])
        if not cands:
            rows.append((insert_id, QUALITY.NO_HITS.name.lower(), False, 0, "", "", ""))
            continue
        scored = []
        for q_val, cid, seq in cands:
            scored.append((_score(q_val, len(seq)), q_val, cid, seq))
        scored.sort(key=lambda t: t[0], reverse=True)
        _, q_val, cid, seq = scored[0]
        paired = q_val in (QUALITY.FULL_MATCH.value, QUALITY.FULL_SCAFFOLD.value, QUALITY.GAPPED_SCAFFOLD.value)
        orig_id = ";".join(id_map.get(c, c) for c in cid.split(";"))
        rows.append((insert_id, _name(q_val), paired, len(seq), "assembly", cid, orig_id))
        out_seq = seq if q_val != QUALITY.REVERSE_END_ONLY.value else str(RawSeq(seq).reverse_complement())
        fa.write(f">{{insert_id}} {{_name(q_val)}} length={{len(seq)}}\\n{{out_seq}}\\n")

cols = ["id", "mapping_quality", "paired", "resolved_length", "assemblers", "contig_id", "original_id"]
pd.DataFrame(rows, columns=cols).to_csv(report_path, index=False)
print(f"resolved={{len(rows)}}")
""")

    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {resolve_script} {ofa.container} {orep.container}",
    )

    return ExecutionResult(
        manifest=[
            {
                out_fa: ofa.local,
                out_rep: orep.local,
            },
        ],
        success=ofa.local.exists() and orep.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=asm,
    resources=Resources(
        cpus=4,
        memory=Size.GB(16),
        duration=Duration(hours=6),
    )
)
