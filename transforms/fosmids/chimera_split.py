"""Split coverage-confirmed two-clone chimeras out of the clustered contigs.

Ports scadc `over_length_and_split.arm_split` + `coverage_profile.per_pool_spans`.
A chimera is two physically distinct clones fused into one contig: a LEFT pool
group covers [0, interior) and a RIGHT pool group covers (interior, L], the two
either overlapping (1x-2x-1x) or abutting at a junction.  Every threshold is
RELATIVE to the contig's own depth -- a pool "covers" a region where its depth
>= split_cover_frac * (that pool's 99th-pctile depth) -- so the caller generalizes
beyond this dataset.  The relative pool-group hard-stop test (each arm carried by a
DISTINCT high-depth pool that drops at the boundary) is what rejected the FF.08B18
false positive in the pilot; keep it.

Interval provenance is baked into each output FASTA header
    >{source}.A source={source} start={s} end={e} action={act}
so the coordinate on the original clustered contig travels with the sequence and
coverage_trim can slice the immutable profile by [start:end].  Uncalled contigs
pass through whole with action=keep.  An optional `verified_chimeras` JSON param
(contig -> {kind, shared|boundary, left_pool, right_pool}) forces specific splits.
"""
from metasmith.python_api import *

lib       = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model     = Transform()
exp       = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
clustered = model.AddRequirement(lib.GetType("fosmids::clustered_contigs"), parents={exp})
profiles  = model.AddRequirement(lib.GetType("fosmids::pool_coverage_profiles"), parents={exp})
img_pyds  = model.AddRequirement(lib.GetType("containers::python_for_data_science.oci"))
out_split = model.AddProduct(lib.GetType("fosmids::split_contigs"))

def protocol(context: ExecutionContext):
    iclustered = context.Input(clustered)
    iprof      = context.Input(profiles)
    osplit     = context.Output(out_split)

    split_cover_frac = context.params.get("split_cover_frac", 0.30)
    split_cover_abs  = context.params.get("split_cover_abs", 5.0)
    min_span_frac    = context.params.get("min_span_frac", 0.15)
    margin           = context.params.get("split_margin", 0.15)
    min_arm_frac     = context.params.get("min_arm_frac", 0.12)
    max_gap_frac     = context.params.get("max_gap_frac", 0.06)
    cover_gap        = context.params.get("cover_gap", 300)
    verified_json    = context.params.get("verified_chimeras", "{}")

    split_script = "chimera_split.py"
    with open(split_script, "w") as f:
        f.write(f"""\
import sys, gzip, pickle, json
import numpy as np
from Bio import SeqIO

FRAC     = {split_cover_frac}
MIN_ABS  = {split_cover_abs}
MIN_SPAN = {min_span_frac}
MARGIN   = {margin}
MIN_ARM  = {min_arm_frac}
MAX_GAP  = {max_gap_frac}
GAP      = {cover_gap}
VERIFIED = json.loads({verified_json!r})

clustered = sys.argv[1]
prof_path = sys.argv[2]
out_path  = sys.argv[3]

with gzip.open(prof_path, "rb") as gz:
    blob = pickle.load(gz)
profiles, lengths = blob["profiles"], blob["lengths"]
seqs = {{s.id: str(s.seq) for s in SeqIO.parse(clustered, "fasta")}}

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

def covered_interval(arr):
    if arr.max() <= 0:
        return None
    peak = np.percentile(arr, 99)
    thr = max(MIN_ABS, FRAC * peak)
    mask = _fill_gaps(arr >= thr, GAP)
    if not mask.any():
        return None
    idx = np.where(mask)[0]
    return int(idx[0]), int(idx[-1] + 1)

def per_pool_spans(contig):
    L = lengths[contig]
    out = {{}}
    for pool, arr in profiles[contig].items():
        iv = covered_interval(arr)
        if iv is None:
            continue
        s, e = iv
        if (e - s) / L < MIN_SPAN:
            continue
        out[pool] = (s, e)
    return out

def arm_split(contig):
    L = lengths[contig]
    spans = per_pool_spans(contig)
    if len(spans) < 2:
        return None
    lo, hi = MARGIN * L, (1 - MARGIN) * L
    left  = [(p, s, e) for p, (s, e) in spans.items() if s <= lo and e < hi]
    right = [(p, s, e) for p, (s, e) in spans.items() if s > lo and e >= hi]
    if not left or not right:
        return None
    lp, ls, le = max(left, key=lambda t: t[2])    # left pool reaching furthest right
    rp, rs, re = min(right, key=lambda t: t[1])   # right pool starting furthest left
    if lp == rp:
        return None                               # arms must be DISTINCT pools
    if le / L < MIN_ARM or (L - rs) / L < MIN_ARM:
        return None
    gap = rs - le                                 # <0 overlap ; >0 junction gap
    if gap > MAX_GAP * L:
        return None
    kind = "overlap" if gap < 0 else "junction"
    boundary = int((le + rs) // 2)
    return dict(kind=kind, left_pool=lp, right_pool=rp,
                left_end=int(le), right_start=int(rs), boundary=boundary)

def segments(contig):
    \"\"\"list of (tag, start, end, action, pool_basis) on the original contig.\"\"\"
    L = lengths[contig]
    if contig in VERIFIED:
        d = VERIFIED[contig]
        if d["kind"] == "overlap":
            os_, oe = d["shared"]
            return [("A", 0, int(oe), "split-overlap", f"{{d.get('left_pool','')}}+shared"),
                    ("B", int(os_), L, "split-overlap", f"{{d.get('right_pool','')}}+shared")]
        bd = int(d["boundary"])
        return [("A", 0, bd, "split-junction", d.get("left_pool", "")),
                ("B", bd, L, "split-junction", d.get("right_pool", ""))]
    d = arm_split(contig)
    if d is None:
        return [(None, 0, L, "keep", "")]
    if d["kind"] == "overlap":
        return [("A", 0, d["left_end"], "split-overlap", f"{{d['left_pool']}}+shared"),
                ("B", d["right_start"], L, "split-overlap", f"{{d['right_pool']}}+shared")]
    bd = d["boundary"]
    return [("A", 0, bd, "split-junction", d["left_pool"]),
            ("B", bd, L, "split-junction", d["right_pool"])]

n_split = 0
with open(out_path, "w") as out:
    for c in sorted(lengths):
        if c not in seqs:
            continue
        segs = segments(c)
        if len(segs) > 1:
            n_split += 1
        for tag, s, e, action, basis in segs:
            piece_id = c if tag is None else f"{{c}}.{{tag}}"
            basis_field = f" pool_basis={{basis}}" if basis else ""
            out.write(f">{{piece_id}} source={{c}} start={{s}} end={{e}} "
                      f"action={{action}}{{basis_field}}\\n{{seqs[c][s:e]}}\\n")

print(f"contigs={{len(lengths)}} chimeras_split={{n_split}}")
""")
    context.ExecWithContainer(
        image=img_pyds,
        cmd=f"python {split_script} {iclustered.container} {iprof.container} {osplit.container}",
    )

    return ExecutionResult(
        manifest=[{out_split: osplit.local}],
        success=osplit.local.exists(),
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
