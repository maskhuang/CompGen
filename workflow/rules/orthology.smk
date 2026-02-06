# Module: orthology
# Rules for orthology analysis (OrthoFinder)
# Implementation: Story 3.1 (OrthoFinder Adapter), Story 3.2 (Orthogroups Table)
#
# This module runs OrthoFinder for orthogroup inference across multiple species,
# then generates structured orthogroup tables for downstream analysis.
# Input: standardized protein sequences from results/standardized/{species}/
# Output: orthogroups, gene counts, species tree, structured tables in results/orthology/

import gzip
import shutil
from pathlib import Path


# =============================================================================
# Helper Functions
# =============================================================================

def get_all_protein_files():
    """Get protein sequence paths for all configured species."""
    return expand(
        "{output_dir}/standardized/{species}/proteins.longest.fa.gz",
        output_dir=get_config_value("project.output_dir", "./results"),
        species=get_species_names(),
    )


# =============================================================================
# Orthology Rules
# =============================================================================

rule orthology_prepare_proteins:
    """
    Decompress protein sequences into a single directory for OrthoFinder input.

    OrthoFinder requires a directory with one .fa file per species (uncompressed).
    File names (without extension) are used as species identifiers.

    AC3: Protein sequence input preparation
    """
    input:
        proteins=get_all_protein_files,
    output:
        proteins_dir=directory("results/orthology/input_proteins"),
        done="results/orthology/input_proteins/.prepared",
    log:
        "logs/orthology_prepare_proteins/run.log",
    run:
        out_dir = Path(output.proteins_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        for proteins_gz in input.proteins:
            # Extract species name from path: .../standardized/{species}/proteins.longest.fa.gz
            species = Path(proteins_gz).parent.name
            out_fa = out_dir / f"{species}.fa"
            with gzip.open(str(proteins_gz), 'rb') as f_in:
                with open(str(out_fa), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        Path(output.done).touch()


rule orthology_infer:
    """
    Run OrthoFinder to infer orthogroups across all species.

    Uses the OrthoFinderAdapter through run_orthofinder.py script.

    AC1: OrthoFinder invocation
    AC2: Output file generation
    AC4: Version compatibility check
    AC5: Error handling
    AC6: Audit record generation
    """
    input:
        proteins_dir="results/orthology/input_proteins",
        done="results/orthology/input_proteins/.prepared",
    output:
        orthogroups="results/orthology/Orthogroups/Orthogroups.tsv",
        gene_count="results/orthology/Orthogroups/Orthogroups.GeneCount.tsv",
        species_tree="results/orthology/Species_Tree/SpeciesTree_rooted.txt",
        run_json="results/meta/orthology_infer/run.run.json",
    params:
        results_dir="results/orthology",
    threads: get_threads("orthology_infer")
    resources:
        mem_mb=65536,
        runtime="12h",
    log:
        "logs/orthology_infer/run.log",
    conda:
        "../envs/orthofinder.yaml"
    script:
        "../scripts/run_orthofinder.py"


rule orthology_parse_orthogroups:
    """
    Parse OrthoFinder output into structured orthogroup tables.

    Generates:
    - orthogroups.tsv: Long-format gene-to-orthogroup mapping
    - orthogroup_stats.tsv: Per-orthogroup statistics
    - species_overlap.tsv: Pairwise species overlap matrix

    Story 3.2: Orthogroups 表生成
    AC1: Long-format table
    AC2: Statistics table
    AC3: Species overlap matrix
    AC5: Audit record
    """
    input:
        orthogroups_raw="results/orthology/Orthogroups/Orthogroups.tsv",
    output:
        orthogroups="results/orthology/orthogroups.tsv",
        stats="results/orthology/orthogroup_stats.tsv",
        overlap="results/orthology/species_overlap.tsv",
        run_json="results/meta/orthology_parse_orthogroups/run.run.json",
    threads: 1
    log:
        "logs/orthology_parse_orthogroups/run.log",
    script:
        "../scripts/build_orthogroup_tables.py"
