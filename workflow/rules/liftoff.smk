# Module: liftoff
# Rules for annotation lifting using Liftoff
# Implementation: Story 5.1, Story 5.2, Story 5.3
#
# This module provides rules for lifting annotations from a reference
# genome to target genomes using the Liftoff tool:
#   results/liftoff/{reference}_to_{target}/
#   ├── lifted_annotation.gff3          # Lifted annotation
#   ├── unmapped_features.txt           # Unmapped features list
#   ├── liftoff_stats.tsv               # Mapping statistics
#   ├── gene_classification.tsv         # Gene classification results (Story 5.2)
#   ├── missing_genes.tsv               # Missing genes list (Story 5.2)
#   └── lifted_annotation_enhanced.gff3 # Enhanced annotation with liftoff attrs (Story 5.3)


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


def get_classification_targets():
    """Get all classification output targets for Story 5.2."""
    targets = []
    for comp in get_liftoff_comparisons():
        ref = comp["reference"]
        target = comp["target"]
        targets.append(f"results/liftoff/{ref}_to_{target}/gene_classification.tsv")
    return targets


def get_enhance_targets():
    """Get all enhanced annotation output targets for Story 5.3."""
    targets = []
    for comp in get_liftoff_comparisons():
        ref = comp["reference"]
        target = comp["target"]
        targets.append(f"results/liftoff/{ref}_to_{target}/lifted_annotation_enhanced.gff3")
    return targets


def get_species_path(species_name: str, file_type: str, default_ext: str) -> str:
    """
    Get file path for a species from config or use default standardized path.

    Args:
        species_name: Species identifier
        file_type: Type of file ('genome' or 'annotation')
        default_ext: Default file extension for standardized path

    Returns:
        Path to the file (may be gzipped or uncompressed)

    Config format:
        liftoff:
          species_paths:
            ncbi_lct:
              genome: /path/to/genome.fna.gz
              annotation: /path/to/annotation.gff.gz
            lab_lct:
              genome: /path/to/assembly.fasta  # uncompressed OK
    """
    species_paths = config.get("liftoff", {}).get("species_paths", {})
    if species_name in species_paths:
        path = species_paths[species_name].get(file_type)
        if path:
            return path
    # Default to standardized path
    if file_type == "annotation":
        return f"results/standardized/{species_name}/annotation.gff3.gz"
    else:
        return f"results/standardized/{species_name}/genome.fa.gz"


def get_reference_gff(wildcards):
    """Get GFF path for reference species (supports custom paths from config)."""
    return get_species_path(wildcards.reference, "annotation", ".gz")


def get_reference_fa(wildcards):
    """Get genome path for reference species (supports custom paths from config)."""
    return get_species_path(wildcards.reference, "genome", ".fa.gz")


def get_target_fa(wildcards):
    """Get genome path for target species (supports custom paths from config)."""
    return get_species_path(wildcards.target, "genome", ".fa.gz")


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
    resources:
        mem_mb=65536,
        runtime="4h",
    log:
        "logs/liftoff_map/{reference}_to_{target}.log"
    envmodules:
        "minimap2/2.26",
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


# =============================================================================
# Gene Classification Rules (Story 5.2)
# =============================================================================

rule liftoff_classify:
    """
    Classify genes based on Liftoff mapping quality.

    Analyzes Liftoff output to classify genes as:
    - present: Successfully mapped with coverage >= threshold
    - uncertain: Partially mapped (coverage < threshold but > 0)
    - missing: Unmapped or zero coverage

    Source: Story 5.2 - 映射质量分析与缺失检测
    """
    input:
        lifted_gff="results/liftoff/{reference}_to_{target}/lifted_annotation.gff3",
        unmapped="results/liftoff/{reference}_to_{target}/unmapped_features.txt",
    output:
        classification="results/liftoff/{reference}_to_{target}/gene_classification.tsv",
        missing="results/liftoff/{reference}_to_{target}/missing_genes.tsv",
        run_json="results/meta/liftoff_classify/reference={reference}_target={target}.run.json",
    params:
        min_coverage=lambda wildcards: config.get("liftoff", {}).get("min_coverage", 0.5),
        min_identity=lambda wildcards: config.get("liftoff", {}).get("min_identity", 0.5),
    threads: 1
    log:
        "logs/liftoff_classify/{reference}_to_{target}.log"
    conda:
        "../envs/liftoff.yaml"
    script:
        "../scripts/classify_absence.py"


rule liftoff_classify_all:
    """
    Run classification for all Liftoff comparisons.

    Triggers classification for all configured reference-target pairs.
    """
    input:
        get_classification_targets()
    output:
        touch("results/liftoff/.classify_complete")


rule liftoff_absence_summary:
    """
    Generate summary of absence detection across all comparisons.

    Aggregates classification statistics from all comparisons into
    a single summary table.

    Source: Story 5.2 AC4 - 批量比较支持
    """
    input:
        classifications=lambda wildcards: [
            f"results/liftoff/{comp['reference']}_to_{comp['target']}/gene_classification.tsv"
            for comp in get_liftoff_comparisons()
        ]
    output:
        summary="results/liftoff/absence_summary.tsv"
    run:
        import csv
        from pathlib import Path
        from workflow.lib.absence_detection import read_classification_stats_from_tsv

        # Summary fieldnames
        ABSENCE_SUMMARY_FIELDS = [
            "reference", "target", "total_genes", "present",
            "uncertain", "missing", "present_rate", "missing_rate"
        ]

        all_summaries = []

        for class_file in input.classifications:
            class_path = Path(class_file)
            if not class_path.exists():
                continue

            # Extract reference and target from path
            # Format: results/liftoff/{reference}_to_{target}/gene_classification.tsv
            dir_name = class_path.parent.name  # e.g., "human_to_mouse_lemur"
            ref, target = dir_name.split("_to_", 1)

            # Use shared function to read and calculate stats
            stats = read_classification_stats_from_tsv(class_path)

            all_summaries.append({
                "reference": ref,
                "target": target,
                **stats,
            })

        # Write summary
        with open(output.summary, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=ABSENCE_SUMMARY_FIELDS, delimiter='\t')
            writer.writeheader()
            if all_summaries:
                writer.writerows(all_summaries)


# =============================================================================
# Annotation Enhancement Rules (Story 5.3)
# =============================================================================

rule liftoff_enhance:
    """
    Enhance lifted annotation with classification attributes.

    Adds liftoff_coverage, liftoff_identity, liftoff_status, source_gene,
    and reference_species attributes to the GFF3. Filters to only include
    genes with 'present' status (or optionally 'uncertain').

    Source: Story 5.3 - 迁移注释输出
    """
    input:
        lifted_gff="results/liftoff/{reference}_to_{target}/lifted_annotation.gff3",
        classification="results/liftoff/{reference}_to_{target}/gene_classification.tsv",
    output:
        enhanced_gff="results/liftoff/{reference}_to_{target}/lifted_annotation_enhanced.gff3",
        run_json="results/meta/liftoff_enhance/reference={reference}_target={target}.run.json",
    params:
        include_uncertain=lambda wildcards: config.get("liftoff", {}).get("include_uncertain", False),
    threads: 1
    log:
        "logs/liftoff_enhance/{reference}_to_{target}.log"
    conda:
        "../envs/liftoff.yaml"
    script:
        "../scripts/enhance_annotation.py"


rule liftoff_enhance_all:
    """
    Run enhancement for all Liftoff comparisons.

    Triggers enhancement for all configured reference-target pairs.

    Source: Story 5.3 AC4 - 批量比较支持
    """
    input:
        get_enhance_targets()
    output:
        touch("results/liftoff/.enhance_complete")
