# lib-transcriptomics

## Archived Intermediates from eguEpdhP

All intermediates are now on Sockeye persistent arc storage at:
**`/arc/project/st-shallam-1/pwy_group/data/porphyridium_purpureum/`**

- `eguEpdhP-intermediates/assembly/` — genome assembly (22MB)
- `eguEpdhP-intermediates/stringtie_gtfs/` — 9 StringTie GTFs (~2.2-2.8MB each)
- `eguEpdhP-intermediates/orfs/` — pprodigal ORFs (7.6MB)
- `eguEpdhP-intermediates/merged_bams/` — 2 merged BAMs (12GB + 20GB, cached from run 62Tq8B53)
- `star_bams/` — 9 STAR BAMs (~4.4GB each)

Local copies also at `~/projects/eguEpdhP-intermediates/` (no BAMs).

## Active Pipeline Run

- **Run ID:** `62Tq8B53` on Sockeye (2026-03-08, fixes: `mode='copy'`, BUSCO bind path)
- **Script:** `tests/manual/run_euk_braker3_pipeline_sockeye.py`
- **Sockeye path:** `/scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53/`
- **Steps:** downloadBuscoLineage, merge_bams, braker3, stringtie_merge, stringtie_quant, gffread_proteins, busco, eggnog_mapper, pydeseq2, stringtie_count_matrix
- **Targets:** gene_count_table, diff_count_table, eggnog_results, busco_results, braker3_gff
- **Fixes applied:** `slurm.nf` publish mode `rellink→copy`; `busco.py` bind path includes `lineages/eukaryota_odb10`
- **Monitor:** `ssh sockeye 'squeue -u txyliu'` and `ssh sockeye 'tail -20 /scratch/st-shallam-1/pwy_group/metasmith/runs/62Tq8B53/_metasmith/logs.latest/nxf.log'`

## Braker3 Container

Switched from `quay.io/biocontainers/braker3:3.0.8` (missing GeneMark) to `teambraker/braker3:v3.0.7.4` which bundles GeneMark-ETP with academic license. `latest` tag fails on Sockeye's apptainer 1.3.1 (OCI manifest bug fixed in 1.3.3+); `v3.0.7.4` uses the older manifest format and pulls successfully. Key env vars in the image:
- `GENEMARK_PATH=/opt/ETP/bin`
- `AUGUSTUS_CONFIG_PATH=/opt/Augustus/config/`
- `AUGUSTUS_BIN_PATH=/opt/Augustus/bin/`
- `AUGUSTUS_SCRIPTS_PATH=/opt/Augustus/scripts/`
