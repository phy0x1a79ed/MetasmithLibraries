"""Trim over-length pieces to their solidly-covered core -> final putative inserts.

Ports scadc `over_length_and_split.classify_over_length` core-extraction.  For each
piece from chimera_split, read its header interval (source, start, end), slice the
IMMUTABLE per-pool profile to profiles[source][:, start:end], sum to a pooled depth
track, and find the contiguous core where pooled depth >= trim_depth_frac * (the
piece's own median depth) -- a RELATIVE threshold, so it generalizes.  A piece is
trimmed to that core only when it is over-length (> HI) and its covered core is
shorter than the piece (the "over-assembly" case); in-window pieces pass whole.

Coordinates stay absolute on the ORIGINAL clustered contig: a trim composes onto the
header's start offset, and the header is rewritten (action += trim).  Finally every
insert is length-banded against the lambda window (LO/HI derived from exp_length /
exp_length_range, default 36500 / 7500 -> 29000-44000).  Nothing is discarded --
short / over44 residuals are written and flagged in the report so the length
distribution is auditable.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
split     = model.AddRequirement(lib.GetType("fosmids::split_contigs"), parents={exp})
profiles  = model.AddRequirement(lib.GetType("fosmids::pool_coverage_profiles"), parents={exp})
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
out_ins   = model.AddProduct(lib.GetType("fosmids::putative_inserts"))
out_rep   = model.AddProduct(lib.GetType("fosmids::putative_insert_report"))

def protocol(context: ExecutionContext):
    isplit = context.Input(split)
    iprof  = context.Input(profiles)
    oins   = context.Output(out_ins)
    orep   = context.Output(out_rep)

    exp_length       = context.params.get("exp_length", 36500)
    exp_length_range = context.params.get("exp_length_range", 7500)
    trim_depth_frac  = context.params.get("trim_depth_frac", 0.20)
    trim_min_abs     = context.params.get("trim_min_abs", 10.0)
    trim_gap         = context.params.get("trim_gap", 300)

    trim_script = "coverage_trim.py"
    with open(trim_script, "w") as f:
        f.write(f"""\
import sys, gzip, pickle
import numpy as np
from Bio import SeqIO

LO = {exp_length} - {exp_length_range}
HI = {exp_length} + {exp_length_range}
FRAC    = {trim_depth_frac}
MIN_ABS = {trim_min_abs}
GAP     = {trim_gap}

split_path = sys.argv[1]
prof_path  = sys.argv[2]
out_fa     = sys.argv[3]
out_rep    = sys.argv[4]

with gzip.open(prof_path, "rb") as gz:
    blob = pickle.load(gz)
profiles, lengths = blob["profiles"], blob["lengths"]

def band(L):
    return "over44" if L > HI else ("in-win" if L >= LO else "short")

def parse_header(desc):
    toks = desc.split()
    piece_id = toks[0]
    kv = {{}}
    for t in toks[1:]:
        if "=" in t:
            k, v = t.split("=", 1)
            kv[k] = v
    return piece_id, kv

def _fill_gaps(mask, max_gap):
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return mask
    m = mask.copy()
    for gi in np.where(np.diff(idx) > 1)[0]:
        a, b = idx[gi], idx[gi + 1]
        if b - a - 1 <= max_gap:
            m[a:b] = True
    return m

def pooled_depth(source, s, e):
    \"\"\"summed per-bp depth over profiles[source][:, s:e] (relative to source array).\"\"\"
    arrs = profiles.get(source, {{}})
    L = e - s
    if L <= 0 or not arrs:
        return np.zeros(max(L, 0), dtype=np.float64)
    acc = np.zeros(lengths[source], dtype=np.float64)
    for arr in arrs.values():
        acc += arr
    return acc[s:e]

def covered_core(depth):
    \"\"\"contiguous [cs, ce) where pooled depth clears the relative floor, or None.\"\"\"
    covered = depth[depth > 0]
    if covered.size == 0:
        return None, 0.0
    base = float(np.median(covered))
    thr = max(FRAC * base, MIN_ABS)
    mask = _fill_gaps(depth >= thr, GAP)
    if not mask.any():
        return None, base
    idx = np.where(mask)[0]
    return (int(idx[0]), int(idx[-1] + 1)), base

rows = []
n_trim = 0
with open(out_fa, "w") as out:
    for rec in SeqIO.parse(split_path, "fasta"):
        piece_id, kv = parse_header(rec.description)
        source = kv.get("source", piece_id)
        s = int(kv.get("start", 0))
        e = int(kv.get("end", len(rec.seq)))
        action = kv.get("action", "keep")
        basis = kv.get("pool_basis", "")
        seq = str(rec.seq)
        L = len(seq)

        depth = pooled_depth(source, s, e)
        core, base = covered_core(depth)

        new_s, new_e, new_seq, new_action = s, e, seq, action
        # low_confidence is reserved for contigs with essentially NO mapped reads
        # (can't assess).  A merely-shallow in-window contig is a real insert and
        # kept; the covered-core floor only gates the over-length TRIM decision.
        if base <= 0:
            new_action = "low_confidence"
        elif L > HI and core is not None:
            cs, ce = core
            if (ce - cs) < L:               # over-assembly: trim to covered core
                new_s, new_e = s + cs, s + ce
                new_seq = seq[cs:ce]
                new_action = action + "+trim" if action != "keep" else "trim"
                n_trim += 1

        nl = len(new_seq)
        out.write(f">{{piece_id}} source={{source}} start={{new_s}} end={{new_e}} "
                  f"action={{new_action}}"
                  + (f" pool_basis={{basis}}" if basis else "")
                  + f"\\n{{new_seq}}\\n")
        rows.append((piece_id, source, new_action, new_s, new_e, nl, basis, band(nl)))

import csv
with open(out_rep, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["insert_id", "source", "action", "start", "end", "length", "pool_basis", "band"])
    w.writerows(rows)

lens = np.array([r[5] for r in rows])
in_win = int(((lens >= LO) & (lens <= HI)).sum())
print(f"inserts={{len(rows)}} trimmed={{n_trim}} "
      f"short={{int((lens < LO).sum())}} in-window={{in_win}} over44={{int((lens > HI).sum())}}")
""")
    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {trim_script} {isplit.container} {iprof.container} {oins.container} {orep.container}",
    )

    return ExecutionResult(
        manifest=[{out_ins: oins.local, out_rep: orep.local}],
        success=oins.local.exists() and orep.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(
        cpus=2,
        memory=Size.GB(16),
        duration=Duration(hours=2),
    )
)
