"""Lane 3 -- DIAMOND blastp of the fosmid ORFs vs UniRef50 -> uniref_hits.

Fresh DIAMOND blastp (mamba, fabfos-annot env) against the staged UniRef50 .dmnd,
emitting BLAST6 + stitle + an appended BSR column (hit_bitscore / self_bitscore),
the 14-col shape `resources/lib/fabfos_evidence.py::read_uniref50` consumes. Self
bitscore is the BLOSUM62 query diagonal converted to bits (DIAMOND KA params) --
no query-vs-query alignment. Ported from scadc functionalAnnotation/diamond_align.py.
The compiler projects the UniProt representative -> MNXR via the Rhea bridge.

GATED: needs the UniRef50 DIAMOND DB (~24GB .dmnd) staged and the fabfos-annot env.
"""
import math
from pathlib import Path
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()
exp   = model.AddRequirement(lib.GetType("fosmids::recovery_experiment"))
orfs  = model.AddRequirement(lib.GetType("sequences::open_reading_frames"), parents={exp})
db    = model.AddRequirement(lib.GetType("functional_annotation::uniref50_dmnd"))
env   = model.AddRequirement(lib.GetType("envs::diamond.condaenv"))
hits  = model.AddProduct(lib.GetType("functional_annotation::uniref_hits"))

_BLOSUM62_DIAG = {
    "A": 4, "R": 5, "N": 6, "D": 6, "C": 9, "Q": 5, "E": 5, "G": 6, "H": 8, "I": 4,
    "L": 4, "K": 5, "M": 5, "F": 6, "P": 7, "S": 4, "T": 5, "W": 11, "Y": 7, "V": 4,
    "B": 4, "Z": 4, "J": 3, "X": 0, "*": 0, "U": 0, "O": 0,
}
_LAMBDA, _K, _LN2 = 0.267, 0.041, math.log(2)
_LN_K = math.log(_K)
_RAW = "diamond_raw.tsv"

def _self_bitscores(fasta_path):
    scores, seq_id, raw = {}, None, 0
    def finalize(sid, r):
        scores[sid] = max(0.0, (_LAMBDA * r - _LN_K) / _LN2)
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if seq_id is not None:
                    finalize(seq_id, raw)
                seq_id, raw = line[1:].split()[0], 0
            elif seq_id is not None:
                for ch in line.strip().upper():
                    raw += _BLOSUM62_DIAG.get(ch, 0)
        if seq_id is not None:
            finalize(seq_id, raw)
    return scores

def protocol(context: ExecutionContext):
    iorfs = context.Input(orfs)
    idb   = context.Input(db)
    ohits = context.Output(hits)
    threads = context.params.get("cpus", 8)
    max_target_seqs = context.params.get("max_target_seqs", 5)
    evalue = context.params.get("evalue", 1e-5)
    context.ExecWithContainer(
        image=env,
        cmd=f"""diamond blastp \
            --query {iorfs.container} \
            --db {idb.container} \
            --out {_RAW} \
            --threads {threads} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
            --max-target-seqs {max_target_seqs} \
            --evalue {evalue} \
            --sensitive""",
    )
    self_bs = _self_bitscores(iorfs.local)
    with open(_RAW) as fin, open(ohits.local, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                bits = float(fields[11])
                sb = self_bs.get(fields[0], 0.0)
                bsr = bits / sb if sb > 0 else 0.0
            except (ValueError, IndexError):
                bsr = 0.0
            fout.write(line.rstrip("\n") + f"\t{bsr:.4f}\n")
    return ExecutionResult(
        manifest=[{hits: ohits.local}],
        success=ohits.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=exp,
    resources=Resources(cpus=8, memory=Size.GB(64), duration=Duration(hours=24)),
)
