# Module: qc
# Rules for quality control (BUSCO annotation completeness assessment)
# Implementation: Story 2.4


def get_proteins_path(wildcards):
    """Get standardized protein sequence path for a species."""
    return f"results/standardized/{wildcards.species}/proteins.longest.fa.gz"


def get_all_busco_summaries():
    """Get all BUSCO summary files for configured species."""
    return expand(
        "results/qc/{species}/busco/short_summary.txt",
        species=get_species_names(),
    )


rule qc_busco:
    """
    Run BUSCO to assess annotation completeness for a single species.

    BUSCO (Benchmarking Universal Single-Copy Orthologs) provides quantitative
    measures for the assessment of genome assembly, gene set, and transcriptome
    completeness based on evolutionarily-informed expectations of gene content.

    Uses BuscoAdapter through run_busco.py script to provide:
    - Version checking (AC4)
    - Error classification (AC5)
    - Audit record generation (AC6)

    Source: Story 2.4 BUSCO QC
    """
    input:
        proteins=get_proteins_path,
    output:
        summary="results/qc/{species}/busco/short_summary.txt",
        full_table="results/qc/{species}/busco/full_table.tsv",
        missing_list="results/qc/{species}/busco/missing_busco_list.tsv",
        run_json="results/meta/qc_busco/species={species}.run.json",
    params:
        lineage=lambda wildcards: get_config_value("busco.lineage", "eukaryota_odb10"),
        mode=lambda wildcards: get_config_value("busco.mode", "proteins"),
        download_path=lambda wildcards: get_config_value("busco.download_path", ""),
        offline=lambda wildcards: get_config_value("busco.offline", False),
        out_dir="results/qc/{species}/busco",
    threads: get_threads("qc_busco")
    log:
        "logs/qc_busco/{species}.log",
    conda:
        "../envs/busco.yaml"
    script:
        "../scripts/run_busco.py"


rule qc_busco_summary:
    """
    Summarize BUSCO results across all species.

    Creates a TSV file with completeness statistics for each species,
    enabling cross-species comparison of annotation quality.

    Source: Story 2.4 BUSCO QC (AC3)
    """
    input:
        summaries=get_all_busco_summaries(),
    output:
        summary="results/qc/busco_summary.tsv",
    log:
        "logs/qc_busco/summary.log",
    script:
        "../scripts/summarize_busco.py"
