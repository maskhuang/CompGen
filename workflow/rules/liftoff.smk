# Module: liftoff
# Rules for annotation lifting using Liftoff
# Implementation: Story 5.1
#
# This module provides rules for lifting annotations from a reference
# genome to target genomes using the Liftoff tool:
#   results/liftoff/{reference}_to_{target}/
#   ├── lifted_annotation.gff3      # Lifted annotation
#   ├── unmapped_features.txt       # Unmapped features list
#   └── liftoff_stats.tsv           # Mapping statistics


# =============================================================================
# Helper Functions
# =============================================================================

def get_liftoff_comparisons():
    """
    Get all Liftoff comparisons from config.

    Returns list of dicts with 'reference' and 'target' keys.

    Config format:
        liftoff:
          comparisons:
            - reference: human
              targets:
                - mouse_lemur
                - ring_tailed_lemur
    """
    comparisons = config.get("liftoff", {}).get("comparisons", [])
    result = []

    for comp in comparisons:
        ref = comp.get("reference")
        targets = comp.get("targets", [])

        if ref and targets:
            for target in targets:
                result.append({"reference": ref, "target": target})

    return result


def get_liftoff_targets():
    """Get all Liftoff output targets for rule all."""
    targets = []
    for comp in get_liftoff_comparisons():
        ref = comp["reference"]
        target = comp["target"]
        targets.append(f"results/liftoff/{ref}_to_{target}/lifted_annotation.gff3")
    return targets


def get_reference_gff(wildcards):
    """Get standardized GFF path for reference species."""
    return f"results/standardized/{wildcards.reference}/annotation.gff3.gz"


def get_reference_fa(wildcards):
    """Get standardized genome path for reference species."""
    return f"results/standardized/{wildcards.reference}/genome.fa.gz"


def get_target_fa(wildcards):
    """Get standardized genome path for target species."""
    return f"results/standardized/{wildcards.target}/genome.fa.gz"


# =============================================================================
# Rules
# =============================================================================

rule liftoff_map:
    """
    Lift annotations from reference genome to target genome using Liftoff.

    Uses LiftoffAdapter through run_liftoff.py script to provide:
    - Version checking (AC4)
    - Error classification (AC5)
    - Audit record generation (AC6)

    Source: Story 5.1 Liftoff Adapter
    """
    input:
        reference_gff=get_reference_gff,
        reference_fa=get_reference_fa,
        target_fa=get_target_fa,
    output:
        lifted_gff="results/liftoff/{reference}_to_{target}/lifted_annotation.gff3",
        unmapped="results/liftoff/{reference}_to_{target}/unmapped_features.txt",
        stats="results/liftoff/{reference}_to_{target}/liftoff_stats.tsv",
        run_json="results/meta/liftoff_map/reference={reference}_target={target}.run.json",
    params:
        min_coverage=lambda wildcards: config.get("liftoff", {}).get("min_coverage", 0.5),
        min_identity=lambda wildcards: config.get("liftoff", {}).get("min_identity", 0.5),
        copies=lambda wildcards: config.get("liftoff", {}).get("copies", False),
        flank=lambda wildcards: config.get("liftoff", {}).get("flank", 0),
        out_dir="results/liftoff/{reference}_to_{target}",
    threads: workflow.cores if hasattr(workflow, 'cores') else 8
    log:
        "logs/liftoff_map/{reference}_to_{target}.log"
    conda:
        "../envs/liftoff.yaml"
    script:
        "../scripts/run_liftoff.py"


rule liftoff_all:
    """
    Run Liftoff for all configured comparisons.

    This is a checkpoint rule that triggers all Liftoff mappings
    defined in the config.
    """
    input:
        get_liftoff_targets()
    output:
        touch("results/liftoff/.liftoff_complete")


rule liftoff_summary:
    """
    Generate summary of all Liftoff runs.

    Aggregates statistics from all Liftoff comparisons into a single TSV file.
    """
    input:
        stats=lambda wildcards: [
            f"results/liftoff/{comp['reference']}_to_{comp['target']}/liftoff_stats.tsv"
            for comp in get_liftoff_comparisons()
        ]
    output:
        summary="results/liftoff/liftoff_summary.tsv"
    run:
        import csv
        from pathlib import Path

        # Define consistent fieldnames for summary output
        SUMMARY_FIELDNAMES = [
            "reference", "target", "lifted_genes", "lifted_features",
            "unmapped_genes", "total_genes", "lift_rate"
        ]

        # Collect all stats
        all_stats = []
        for stats_file in input.stats:
            stats_path = Path(stats_file)
            if stats_path.exists():
                with open(stats_path, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        all_stats.append(row)

        # Write summary with consistent fieldnames
        with open(output.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDNAMES, delimiter='\t')
            writer.writeheader()
            if all_stats:
                writer.writerows(all_stats)
