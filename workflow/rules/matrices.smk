# Module: matrices
# Rules for presence/absence and copy number matrix generation
# Implementation: Story 6A.1 (Presence/Absence Matrix), Story 6A.2 (Multi-species Comparison)
#
# This module generates matrices from orthogroup data for downstream analysis.
# Input: structured orthogroups table from results/orthology/
# Output: matrices in results/matrices/


# =============================================================================
# Matrix Generation Rules
# =============================================================================

rule matrices_generate:
    """
    Generate presence/absence and copy number matrices from orthogroups.

    Reads the long-format orthogroups.tsv and produces:
    - presence_absence.tsv: Binary 0/1 matrix (rows=orthogroups, cols=species)
    - copy_number.tsv: Gene count matrix (rows=orthogroups, cols=species)

    Story 6A.1: FR20
    AC1: Binary presence/absence matrix
    AC2: Copy number matrix
    AC3: Matrix consistency
    AC5: Audit record
    """
    input:
        orthogroups="results/orthology/orthogroups.tsv",
    output:
        presence_absence="results/matrices/presence_absence.tsv",
        copy_number="results/matrices/copy_number.tsv",
        run_json="results/meta/matrices_generate/run.run.json",
    threads: 1
    log:
        "logs/matrices_generate/run.log",
    script:
        "../scripts/build_matrices.py"
