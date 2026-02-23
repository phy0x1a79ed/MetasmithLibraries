# Plan: Update Binning Workflow Tests to Use AB48 Community Data

## Context
The binning workflow tests currently reference `small_assembly.fna` and `small_bam.bam` which don't exist. We need to use the real AB48 community data (`ABC-240403_KD.fna` assembly + `ABC-240403_KD.fastq.gz` reads) and generate a BAM file by running the Metasmith `assembly_stats` transform via Docker.

## Files to Modify
- `tests/conftest.py` — add BAM-generation fixture and AB48-specific fixtures
- `tests/test_binning_workflow.py` — update `binning_input` to use AB48 data

## Approach

### 1. Add AB48 BAM generation fixture to `conftest.py`

Add a **session-scoped** fixture `ab48_bam` that:
1. Checks if `tests/test_data/ab48_community/ABC-240403_KD.bam` already exists → returns it immediately (cached)
2. If not, runs the `assembly_stats` workflow via the Metasmith agent to produce the BAM:
   - Creates an input library with:
     - `sequences::read_metadata` (AddValue: `{"parity": "single", "length_class": "long"}`)
     - `sequences::long_reads` (the fastq.gz, parents={meta})
     - `sequences::assembly` (the .fna, parents={reads})
   - Loads assembly transforms from `transforms/assembly`
   - Targets `alignment::bam` (the solver will auto-include read_qc_stats generation as an intermediate step since assembly_stats requires it)
   - Generates, stages, and runs the workflow
   - Copies the resulting BAM to `tests/test_data/ab48_community/ABC-240403_KD.bam`

This is a heavy operation (~1.6GB reads) so caching the BAM to disk is essential.

### 2. Update `binning_input` fixture in `test_binning_workflow.py`

Change `binning_input` to:
- Reference `test_data/ab48_community/ABC-240403_KD.fna` as the assembly
- Depend on the `ab48_bam` fixture for the BAM file
- Keep the same AddItem pattern: assembly (no parents), BAM (parents={asm})

### 3. Key details
- Use `sequences::long_reads` for the reads
- `read_metadata`: `{"parity": "single", "length_class": "long"}`
- The solver will plan read_qc_stats generation automatically since assembly_stats requires it
- Assembly transforms loaded from `mlib / "transforms/assembly"`
- Resource overrides: 14 CPUs, enough memory for minimap2 on 1.6GB reads

## Verification
1. Run `pytest tests/test_binning_workflow.py::TestBinningWorkflowGeneration::test_can_plan_metabat2_workflow -v` to verify planning works
2. Run the full binning E2E tests with `pytest tests/test_binning_workflow.py -v --slow` to verify end-to-end
