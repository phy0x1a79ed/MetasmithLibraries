"""
Planning-path coverage matrix across transform libraries.

These tests validate that representative workflow paths can be generated
for all transform libraries without executing workflows.
"""
from dataclasses import dataclass
from pathlib import Path

import pytest
from metasmith.python_api import (
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
)


@dataclass(frozen=True)
class PathCase:
    id: str
    transforms: tuple[str, ...]
    target: str
    sample_type: str
    profile: str
    needs_lib_resources: bool = False
    allow_zero_steps: bool = False


CASES = [
    # assembly
    PathCase("assembly-read-qc", ("assembly",), "sequences::read_qc_stats", "sequences::read_metadata", "assembly_short"),
    PathCase("assembly-megahit", ("assembly",), "sequences::assembly", "sequences::read_metadata", "assembly_short"),
    PathCase("assembly-flye", ("assembly",), "sequences::flye_assembly", "sequences::long_reads", "assembly_long"),
    PathCase("assembly-hifiasm-meta", ("assembly",), "sequences::hifiasm_meta_assembly", "sequences::long_reads", "assembly_long"),
    # amplicon
    PathCase("amplicon-asv-map", ("amplicon",), "amplicon::asv_contig_map", "amplicon::asv_seqs", "amplicon"),
    PathCase("amplicon-asv-taxonomy", ("amplicon",), "amplicon::asv_taxonomy", "amplicon::asv_seqs", "amplicon"),
    # functional annotation
    PathCase("annotation-pprodigal", ("functionalAnnotation",), "sequences::orfs", "sequences::assembly", "annotation_assembly"),
    PathCase("annotation-interproscan", ("functionalAnnotation",), "annotation::interproscan_json", "sequences::orfs", "annotation_orfs"),
    PathCase("annotation-kofamscan", ("functionalAnnotation",), "annotation::kofamscan_results", "sequences::orfs", "annotation_orfs"),
    PathCase("annotation-proteinbert", ("functionalAnnotation",), "annotation::proteinbert_embeddings", "sequences::orfs", "annotation_orfs"),
    PathCase("annotation-eggnog", ("functionalAnnotation",), "annotation::eggnog_mapper_results", "sequences::orfs", "annotation_orfs"),
    PathCase("annotation-busco", ("functionalAnnotation",), "annotation::busco_results", "transcriptomics::braker3_proteins", "annotation_busco"),
    PathCase("annotation-virsorter2", ("functionalAnnotation",), "annotation::virsorter2_viral_sequences", "sequences::assembly", "annotation_assembly"),
    PathCase("annotation-dramv", ("functionalAnnotation",), "annotation::dramv_annotations", "sequences::assembly", "annotation_assembly"),
    # logistics
    PathCase("logistics-get-assembly", ("logistics",), "sequences::assembly", "ncbi::assembly_accession", "logistics_ncbi"),
    PathCase("logistics-get-sra-cache", ("logistics",), "ncbi::sra_cache", "sequences::read_metadata", "logistics_ncbi"),
    PathCase("logistics-dump-sra-short", ("logistics",), "sequences::short_reads", "sequences::read_metadata", "logistics_ncbi"),
    PathCase("logistics-download-kofam", ("logistics",), "annotation::kofamscan_profiles", "annotation::kofamscan_source", "logistics_sources"),
    PathCase("logistics-download-silva", ("logistics",), "amplicon::silva_db", "amplicon::silva_source", "logistics_sources"),
    PathCase("logistics-interleave-short", ("logistics",), "sequences::short_reads", "sequences::read_pair", "logistics_pair_plain"),
    PathCase("logistics-ora2fastq", ("logistics",), "sequences::short_reads", "sequences::read_pair", "logistics_pair_ora"),
    PathCase("logistics-pull-container", ("logistics",), "containers::pulled_container", "containers::container", "logistics_container"),
    # metagenomics
    PathCase("metagenomics-checkm", ("metagenomics",), "taxonomy::checkm_stats", "sequences::assembly", "metagenomics_taxonomy"),
    PathCase("metagenomics-gtdbtk", ("metagenomics",), "taxonomy::gtdbtk", "sequences::assembly", "metagenomics_taxonomy"),
    PathCase("metagenomics-metabuli", ("metagenomics",), "taxonomy::metabuli", "sequences::assembly", "metagenomics_taxonomy"),
    PathCase("metagenomics-fastani", ("metagenomics",), "taxonomy::ani_table", "pangenome::pangenome", "metagenomics_fastani"),
    PathCase("metagenomics-metabat2", ("metagenomics",), "binning::metabat2_bin_fasta", "sequences::assembly", "metagenomics_binning"),
    PathCase("metagenomics-semibin2", ("metagenomics",), "binning::semibin2_bin_fasta", "sequences::assembly", "metagenomics_binning"),
    PathCase("metagenomics-comebin", ("metagenomics",), "binning::comebin_bin_fasta", "sequences::assembly", "metagenomics_binning"),
    # pangenome
    PathCase("pangenome-matrix", ("pangenome",), "pangenome::ppanggolin_matrix", "pangenome::pangenome", "pangenome"),
    PathCase("pangenome-heatmap", ("pangenome",), "pangenome::heatmap", "pangenome::pangenome", "pangenome", needs_lib_resources=True),
    # response surface
    PathCase(
        "response-surface-coefficients",
        ("responseSurface",),
        "media_optimization::response_surface_coefficients",
        "media_optimization::growth_data",
        "response_surface",
        needs_lib_resources=True,
    ),
    # transcriptomics
    PathCase("transcriptomics-salmon-quant", ("transcriptomics",), "transcriptomics::salmon_quant", "transcriptomics::experiment", "transcriptomics"),
    PathCase("transcriptomics-count-table", ("transcriptomics",), "transcriptomics::count_table", "transcriptomics::experiment", "transcriptomics"),
    PathCase("transcriptomics-gene-count", ("transcriptomics",), "transcriptomics::gene_count_table", "transcriptomics::experiment", "transcriptomics"),
    PathCase("transcriptomics-diff-count", ("transcriptomics",), "transcriptomics::diff_count_table", "transcriptomics::experiment", "transcriptomics"),
    PathCase("transcriptomics-volcano", ("transcriptomics",), "transcriptomics::volcano_plot", "transcriptomics::experiment", "transcriptomics"),
    # metabolic modelling
    PathCase(
        "metabolic-feature-table",
        ("metabolicModelling",),
        "metabolomics::metabolomics_feature_table",
        "metabolomics::jgi_metabolomics_dataset",
        "metabolomics",
        allow_zero_steps=True,
    ),
    PathCase(
        "metabolic-differential",
        ("metabolicModelling",),
        "metabolomics::metabolomics_differential",
        "metabolomics::jgi_metabolomics_dataset",
        "metabolomics",
        allow_zero_steps=True,
    ),
]


@pytest.fixture(scope="session")
def load_transform_library(mlib):
    cache = {}

    def _load(name: str):
        if name not in cache:
            cache[name] = TransformInstanceLibrary.Load(mlib / "transforms" / name)
        return cache[name]

    return _load


def _ensure_text(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(text)
    return path


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def _new_inputs(tmp_path: Path, mlib: Path, profile: str) -> DataInstanceLibrary:
    inputs = DataInstanceLibrary(tmp_path / f"{profile}.xgdb")
    for type_lib in sorted((mlib / "data_types").glob("*.yml")):
        inputs.AddTypeLibrary(type_lib)
    return inputs


def _build_profile_inputs(profile: str, tmp_path: Path, mlib: Path, test_data_dir: Path) -> DataInstanceLibrary:
    scratch = tmp_path / "matrix_data"
    scratch.mkdir(parents=True, exist_ok=True)
    inputs = _new_inputs(tmp_path, mlib, profile)

    small_assembly = _first_existing([test_data_dir / "small_assembly.fna"])
    if small_assembly is None:
        small_assembly = _ensure_text(scratch / "small_assembly.fna", ">ctg1\nATGC\n")

    small_orfs = _first_existing([test_data_dir / "small_orfs.faa"])
    if small_orfs is None:
        small_orfs = _ensure_text(scratch / "small_orfs.faa", ">orf1\nMKKL\n")

    r1_zipped = _first_existing(sorted((test_data_dir / "porphyridium" / "Raw-transcriptome-data").glob("*_R1_001.fastq.gz")))
    r2_zipped = _first_existing(sorted((test_data_dir / "porphyridium" / "Raw-transcriptome-data").glob("*_R2_001.fastq.gz")))
    if r1_zipped is None:
        r1_zipped = _ensure_text(scratch / "reads_R1.fastq.gz", "@r1\nACGT\n+\n####\n")
    if r2_zipped is None:
        r2_zipped = _ensure_text(scratch / "reads_R2.fastq.gz", "@r2\nTGCA\n+\n####\n")

    long_reads = _first_existing([
        test_data_dir / "ab48_community" / "ABC-240403_KD.fastq.gz",
        test_data_dir / "mock_reads.fq",
    ])
    if long_reads is None:
        long_reads = _ensure_text(scratch / "long_reads.fq.gz", "@lr\nACGTACGT\n+\n########\n")

    small_bam = _first_existing([
        test_data_dir / "ab48_community" / "ABC-240403_KD.bam",
        test_data_dir / "small_bam.bam",
    ])
    if small_bam is None:
        small_bam = _ensure_text(scratch / "small.bam", "BAM\n")

    kofam_profiles = _first_existing([test_data_dir / "kofam" / "profiles"])
    if kofam_profiles is None:
        kofam_profiles = _ensure_dir(scratch / "kofam_profiles")
        _ensure_text(kofam_profiles / "K00001.hmm", "HMMER3\n")
    kofam_ko_list = _first_existing([test_data_dir / "kofam" / "ko_list"])
    if kofam_ko_list is None:
        kofam_ko_list = _ensure_text(scratch / "ko_list.tsv", "K00001\tmock\n")

    interpro_data = _first_existing([test_data_dir / "interproscan_data"])
    if interpro_data is None:
        interpro_data = _ensure_dir(scratch / "interproscan_data")

    uniref_dmnd = _first_existing([test_data_dir / "uniref50" / "uniref50.dmnd"])
    if uniref_dmnd is None:
        uniref_dmnd = _ensure_text(scratch / "uniref50.dmnd", "DMND\n")

    ora_r1 = _first_existing([test_data_dir / "ora2fastq" / "S485_S64_L006_R1_001.fastq.ora"])
    ora_r2 = _first_existing([test_data_dir / "ora2fastq" / "S485_S64_L006_R2_001.fastq.ora"])
    if ora_r1 is None:
        ora_r1 = _ensure_text(scratch / "reads_R1.fastq.ora", "ora\n")
    if ora_r2 is None:
        ora_r2 = _ensure_text(scratch / "reads_R2.fastq.ora", "ora\n")

    growth_csv = _first_existing(sorted((test_data_dir / "media_opt" / "data" / "dse").glob("*.csv")))
    if growth_csv is None:
        growth_csv = _ensure_text(
            scratch / "growth.csv",
            "time,replicate,od\n0,A,0.1\n1,A,0.2\n",
        )

    jgi_dataset = _first_existing([test_data_dir / "porphyridium-metabolomics" / "jgi_metabolomic_data"])
    if jgi_dataset is None:
        jgi_dataset = _ensure_dir(scratch / "jgi_metabolomic_data")
        _ensure_dir(jgi_dataset / "Targeted")
        _ensure_dir(jgi_dataset / "Untargeted")

    if profile == "assembly_short":
        meta = inputs.AddValue(
            "reads_metadata_short.json",
            {"parity": "single", "length_class": "short"},
            "sequences::read_metadata",
        )
        short_reads = _ensure_text(scratch / "short_reads.fq.gz", "@r\nACGT\n+\n####\n")
        inputs.AddItem(short_reads, "sequences::short_reads_se", parents={meta})

    elif profile == "assembly_long":
        meta = inputs.AddValue(
            "reads_metadata_long.json",
            {"parity": "single", "length_class": "long"},
            "sequences::read_metadata",
        )
        inputs.AddItem(long_reads, "sequences::long_reads", parents={meta})

    elif profile == "amplicon":
        asvs = _ensure_text(scratch / "small_asvs.fasta", ">asv1\nATGC\n")
        inputs.AddItem(asvs, "amplicon::asv_seqs")
        inputs.AddItem(small_assembly, "sequences::assembly")
        inputs.AddValue("blast_identity_threshold.txt", "97", "amplicon::blast_identity_threshold")
        classifier = _ensure_text(scratch / "silva_classifier.qza", "classifier\n")
        inputs.AddItem(classifier, "amplicon::silva_nb_classifier")

    elif profile == "annotation_orfs":
        inputs.AddItem(small_orfs, "sequences::orfs")
        inputs.AddItem(interpro_data, "annotation::interproscan_data")
        inputs.AddItem(kofam_profiles, "annotation::kofamscan_profiles")
        inputs.AddItem(kofam_ko_list, "annotation::kofamscan_ko_list")
        inputs.AddItem(uniref_dmnd, "annotation::uniref50_diamond_db")
        inputs.AddItem(_ensure_dir(scratch / "eggnog_db"), "annotation::eggnog_db")

    elif profile == "annotation_assembly":
        inputs.AddItem(small_assembly, "sequences::assembly")
        inputs.AddItem(_ensure_dir(scratch / "virsorter2_db"), "annotation::virsorter2_db")
        inputs.AddItem(_ensure_dir(scratch / "dram_db"), "annotation::dram_db")

    elif profile == "annotation_busco":
        proteins = _ensure_text(scratch / "braker3_proteins.faa", ">p1\nMKKL\n")
        inputs.AddItem(proteins, "transcriptomics::braker3_proteins")
        inputs.AddItem(_ensure_dir(scratch / "busco_lineage"), "annotation::busco_lineage")

    elif profile == "logistics_ncbi":
        meta = inputs.AddValue(
            "reads_metadata.json",
            {"parity": "paired", "length_class": "short"},
            "sequences::read_metadata",
        )
        inputs.AddValue("sra_accession.txt", "SRR5585544", "ncbi::sra_accession", parents={meta})
        inputs.AddValue("assembly_accession.txt", "GCF_000005845.2", "ncbi::assembly_accession")
        inputs.AddValue("accession_list.csv", "GCF_000005845.2\nGCA_008690995.1\n", "ncbi::accession_list")

    elif profile == "logistics_sources":
        inputs.AddValue("kofam_source.txt", "ftp://example/kofam", "annotation::kofamscan_source")
        inputs.AddValue("silva_source.txt", "SILVA_138.2_SSURef_NR99", "amplicon::silva_source")

    elif profile == "logistics_pair_plain":
        pair = inputs.AddValue("read_pair_plain.txt", "pair-plain", "sequences::read_pair")
        r1 = _ensure_text(scratch / "plain_R1.fq", "@r1\nACGT\n+\n####\n")
        r2 = _ensure_text(scratch / "plain_R2.fq", "@r2\nTGCA\n+\n####\n")
        inputs.AddItem(r1, "sequences::forward_short_reads", parents={pair})
        inputs.AddItem(r2, "sequences::reverse_short_reads", parents={pair})

    elif profile == "logistics_pair_zipped":
        pair = inputs.AddValue("read_pair_zip.txt", "pair-zipped", "sequences::read_pair")
        inputs.AddItem(r1_zipped, "sequences::zipped_forward_short_reads", parents={pair})
        inputs.AddItem(r2_zipped, "sequences::zipped_reverse_short_reads", parents={pair})

    elif profile == "logistics_pair_ora":
        pair = inputs.AddValue("read_pair_ora.txt", "pair-ora", "sequences::read_pair")
        inputs.AddItem(ora_r1, "sequences::forward_ora_reads", parents={pair})
        inputs.AddItem(ora_r2, "sequences::reverse_ora_reads", parents={pair})

    elif profile == "logistics_container":
        container_file = _ensure_text(scratch / "dummy.oci", "oci\n")
        inputs.AddItem(container_file, "containers::container")

    elif profile == "metagenomics_taxonomy":
        inputs.AddItem(small_assembly, "sequences::assembly")
        inputs.AddItem(_ensure_dir(scratch / "gtdb"), "taxonomy::gtdb")
        inputs.AddItem(_ensure_dir(scratch / "metabuli_ref"), "taxonomy::metabuli_ref")

    elif profile == "metagenomics_fastani":
        pan = inputs.AddValue("pan.txt", "pan-1", "pangenome::pangenome")
        inputs.AddItem(small_assembly, "sequences::assembly", parents={pan})

    elif profile == "metagenomics_binning":
        asm = inputs.AddItem(small_assembly, "sequences::assembly")
        inputs.AddItem(small_bam, "alignment::bam", parents={asm})

    elif profile == "pangenome":
        pan = inputs.AddValue("pangenome.txt", "pan-1", "pangenome::pangenome")
        gbk1 = _ensure_text(scratch / "g1.gbk", "LOCUS       g1\nORIGIN\n//\n")
        gbk2 = _ensure_text(scratch / "g2.gbk", "LOCUS       g2\nORIGIN\n//\n")
        inputs.AddItem(gbk1, "sequences::gbk", parents={pan})
        inputs.AddItem(gbk2, "sequences::gbk", parents={pan})

    elif profile == "response_surface":
        inputs.AddItem(growth_csv, "media_optimization::growth_data")

    elif profile == "transcriptomics":
        exp = inputs.AddValue("experiment.txt", "exp-1", "transcriptomics::experiment")
        pair = inputs.AddValue("pair.txt", "pair-1", "sequences::read_pair", parents={exp})
        inputs.AddItem(r1_zipped, "sequences::zipped_forward_short_reads", parents={pair})
        inputs.AddItem(r2_zipped, "sequences::zipped_reverse_short_reads", parents={pair})
        inputs.AddItem(small_assembly, "sequences::assembly", parents={exp})

    elif profile == "metabolomics":
        inputs.AddItem(jgi_dataset, "metabolomics::jgi_metabolomics_dataset")

    else:
        raise ValueError(f"Unknown matrix profile: {profile}")

    inputs.Save()
    return inputs


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_can_plan_transform_paths_matrix(
    case,
    agent,
    base_resources,
    lib_resources,
    load_transform_library,
    mlib,
    test_data_dir,
    tmp_path,
):
    if case.needs_lib_resources and lib_resources is None:
        pytest.skip(f"{case.id} requires resources/lib")

    inputs = _build_profile_inputs(case.profile, tmp_path, mlib, test_data_dir)
    samples = list(inputs.AsSamples(case.sample_type))
    if not samples:
        pytest.skip(f"No samples available for type {case.sample_type}")

    resources = list(base_resources) + [inputs]
    if lib_resources is not None:
        resources.append(lib_resources)

    targets = TargetBuilder()
    targets.Add(case.target)

    transforms = [load_transform_library(name) for name in case.transforms]
    task = agent.GenerateWorkflow(
        samples=samples,
        resources=resources,
        transforms=transforms,
        targets=targets,
    )

    assert task.ok, f"{case.id}: workflow generation failed for target {case.target}: {task}"
    if not case.allow_zero_steps:
        assert len(task.plan.steps) > 0, f"{case.id}: generated workflow should include at least one step"
