# Current State & Next Steps

## What was done

### 1. Fixed result persistence (slurm.nf)
- **File:** `/home/tony/workspace/msm/Metasmith/src/metasmith/nextflow_config/slurm.nf` line 75
- Changed `mode = 'rellink'` → `mode = 'copy'` so results survive `cleanup = true`

### 2. Fixed BUSCO bind path (busco.py)
- **File:** `transforms/functionalAnnotation/busco.py` line 19
- Changed bind target from `/busco_lineage` to `/busco_lineage/lineages/eukaryota_odb10`
- BUSCO expects `<download_path>/lineages/<lineage>/`, so the lineage dir must be mounted at that nested path

### 3. Cached merged BAMs to arc
- Copied from run `62Tq8B53` work dirs to `eguEpdhP-intermediates/merged_bams/`
- Later re-merged all 9 into single BAM at `/scratch/st-shallam-1/pwy_group/metasmith/cache/merged_bams/porphyridium_all9_iNlpm1XR.bam` (32GB)

### 4. Run mHAnaWSi — braker3 + stringtie_assemble succeeded, pipeline stuck
- **6-step DAG:** stringtie_assemble → braker3 → stringtie_merge → stringtie_quant → pydeseq2 + stringtie_count_matrix
- stringtie_assemble (9x): all COMPLETED exit 0 (3-6 min each)
- braker3 (1x): COMPLETED exit 0 (1h 20m)
- stringtie_merge: **NEVER SUBMITTED** — orchestrator channel deadlock
- Results downloaded to `./data/mHAnaWSi-results/` (9 GTFs + GFF3 + proteins)
- Logs at `./data/mHAnaWSi-results/logs/`

## MetaSmith Bugs

### Bug 1: Batching splits group_by groups (run 62Tq8B53)
MetaSmith split 9 BAMs into two merge_bams groups (6 + 3) instead of one. `group_by=exp` should have grouped all 9 under the single experiment. Logs: `./data/metasmith_batch_bug/`

### Bug 2: Orchestrator channel deadlock (run mHAnaWSi)
After braker3 completes, its output (`_mR9ozKye` = braker3_gff) is not forwarded to stringtie_merge's input group channel. Nextflow reports "No more task to compute" while p03-p06 remain ACTIVE with open ports. The workflow.nf wiring looks correct:
```
(_mR9ozKye, _mnweZ2NO) = o.post([*p02__braker3(o.group('Rsoo8XE7', [...], k, 1))], k)
(_jDFH0gS6) = o.post([*p03__stringtie_merge(o.group('Rsoo8XE7', [_Rsoo8XE7, _vyhx3Bb8, _mR9ozKye, _rRC4FxHn], k, 1))], k)
```
The `o.group` for stringtie_merge groups by `Rsoo8XE7` (experiment) and waits for `_vyhx3Bb8` (9 stringtie GTFs) + `_mR9ozKye` (1 braker3 GFF). Suspect the Orchestrator's channel join logic deadlocks when the cardinalities differ (9 GTFs vs 1 GFF vs 1 experiment).

**Reproduces on resume:** Restarted via `start.sh` — Nextflow cached all 10 tasks instantly but deadlocked again at the same point. The bug is deterministic and not a race condition. Logs: `./data/mHAnaWSi-results/logs/`

**Graceful stop:** delete `PID.lock` in the run directory to stop MetaSmith/Nextflow cleanly.

## What to do next

1. **Fix orchestrator channel deadlock** in MetaSmith — this is blocking all downstream steps from braker3
2. **Workaround:** run stringtie_merge → stringtie_quant → pydeseq2 + count_matrix manually using the downloaded braker3 GFF + stringtie GTFs
3. **Then:** add functional annotation (busco, eggnog_mapper) once count tables work

## Available intermediates (local)
- `./data/mHAnaWSi-results/*.gff3` — braker3 GFF (4.5MB)
- `./data/mHAnaWSi-results/*.faa` — braker3 proteins (4.7MB)
- `./data/mHAnaWSi-results/*.gtf` — 9 stringtie GTFs (318-437KB each)

## Key paths
- **Sockeye scratch:** `/scratch/st-shallam-1/pwy_group/metasmith/runs/`
- **Arc intermediates:** `/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/eguEpdhP-intermediates/`
- **Merged BAM cache:** `/scratch/st-shallam-1/pwy_group/metasmith/cache/merged_bams/`
- **SLURM account secret:** `secrets/slurm_account_sockeye`
- **Monitor:** `ssh sockeye 'squeue -u txyliu'`
