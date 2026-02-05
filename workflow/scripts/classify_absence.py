"""
Classify genes based on Liftoff mapping quality.

This script bridges Snakemake with the absence detection module, providing:
- Gene classification (present/uncertain/missing) based on coverage/identity
- Missing genes list generation
- Audit record generation

Source: Story 5.2 - 映射质量分析与缺失检测

Usage (via Snakemake):
    snakemake.input.lifted_gff: Path to Liftoff lifted_annotation.gff3
    snakemake.input.unmapped: Path to Liftoff unmapped_features.txt
    snakemake.output.classification: Path to gene_classification.tsv
    snakemake.output.missing: Path to missing_genes.tsv
    snakemake.output.run_json: Path to .run.json audit file
    snakemake.params.min_coverage: Minimum coverage threshold
    snakemake.params.min_identity: Minimum identity threshold
    snakemake.wildcards.reference: Reference species name
    snakemake.wildcards.target: Target species name
"""

import csv
import shutil
import time
from pathlib import Path

from workflow.lib.absence_detection import (
    __version__ as absence_detection_version,
    classify_liftoff_genes,
    calculate_classification_stats,
    filter_missing_genes,
)
from workflow.lib.audit import create_and_write_audit
from workflow.lib.errors import CompGeneError


# =============================================================================
# Output Writers
# =============================================================================

# Output field definitions
CLASSIFICATION_FIELDS = [
    "gene_id", "gene_name", "reference_species", "target_species",
    "status", "coverage", "identity", "source_gene"
]

MISSING_FIELDS = [
    "gene_id", "gene_name", "reference_species", "target_species",
    "liftoff_status", "coverage", "identity"
]


def write_classification_tsv(results: list[dict], output_path: Path) -> None:
    """
    Write gene classification results to TSV file.

    Args:
        results: List of classification result dictionaries.
        output_path: Path to write TSV file.
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CLASSIFICATION_FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)


def write_missing_tsv(results: list[dict], output_path: Path) -> None:
    """
    Write missing genes list to TSV file.

    Args:
        results: List of missing gene dictionaries.
        output_path: Path to write TSV file.
    """
    # Transform to missing genes format (rename 'status' to 'liftoff_status')
    missing_rows = []
    for r in results:
        missing_rows.append({
            "gene_id": r["gene_id"],
            "gene_name": r["gene_name"],
            "reference_species": r["reference_species"],
            "target_species": r["target_species"],
            "liftoff_status": "unmapped" if r["status"] == "missing" else "partial",
            "coverage": r["coverage"],
            "identity": r["identity"],
        })

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=MISSING_FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(missing_rows)


# =============================================================================
# Main Execution
# =============================================================================

def run_classification():
    """
    Execute gene classification from Liftoff outputs.

    This function:
    1. Reads Liftoff lifted GFF and unmapped features
    2. Classifies genes based on coverage/identity thresholds
    3. Generates gene_classification.tsv (all genes)
    4. Generates missing_genes.tsv (missing + uncertain genes)
    5. Generates .run.json audit record
    """
    start_time = time.time()

    # Get Snakemake variables
    lifted_gff_path = Path(snakemake.input.lifted_gff)
    unmapped_path = Path(snakemake.input.unmapped)

    classification_output = Path(snakemake.output.classification)
    missing_output = Path(snakemake.output.missing)
    run_json_output = Path(snakemake.output.run_json)

    reference = snakemake.wildcards.reference
    target = snakemake.wildcards.target

    min_coverage = snakemake.params.min_coverage
    min_identity = snakemake.params.min_identity

    # Track execution for audit
    exit_code = 0
    error_code_str = None
    error_message = None
    stats = {}

    try:
        # Step 1: Classify genes
        print(f"Classifying genes from {lifted_gff_path}")
        print(f"Thresholds: min_coverage={min_coverage}, min_identity={min_identity}")

        results = classify_liftoff_genes(
            lifted_gff=lifted_gff_path,
            unmapped_file=unmapped_path,
            reference_species=reference,
            target_species=target,
            min_coverage=min_coverage,
            min_identity=min_identity,
        )

        # Step 2: Calculate statistics
        stats = calculate_classification_stats(results)
        print(f"Classification results: "
              f"{stats['present']} present, "
              f"{stats['uncertain']} uncertain, "
              f"{stats['missing']} missing "
              f"(total: {stats['total_genes']})")

        # Step 3: Write classification output
        classification_output.parent.mkdir(parents=True, exist_ok=True)
        write_classification_tsv(results, classification_output)
        print(f"Gene classification written to: {classification_output}")

        # Step 4: Write missing genes output (includes uncertain)
        missing_genes = filter_missing_genes(results, include_uncertain=True)
        write_missing_tsv(missing_genes, missing_output)
        print(f"Missing genes list written to: {missing_output}")
        print(f"Missing genes count: {len(missing_genes)} "
              f"(missing={stats['missing']}, uncertain={stats['uncertain']})")

    except CompGeneError as e:
        error_code_str = e.error_code.value
        error_message = str(e)
        exit_code = e.to_exit_code()
        raise

    except Exception as e:
        error_code_str = "E_RUNTIME"
        error_message = str(e)
        exit_code = 1
        raise

    finally:
        # Step 5: Write audit record
        runtime = time.time() - start_time

        # Ensure meta dir exists
        run_json_output.parent.mkdir(parents=True, exist_ok=True)

        # Build command with parameters for audit traceability
        cmd = [
            "classify_absence.py",
            f"--min-coverage={min_coverage}",
            f"--min-identity={min_identity}",
            f"--reference={reference}",
            f"--target={target}",
        ]

        _, audit_path = create_and_write_audit(
            rule="liftoff_classify",
            wildcards={"reference": reference, "target": target},
            cmd=cmd,
            tool_version=absence_detection_version,
            input_paths={
                "lifted_gff": lifted_gff_path,
                "unmapped": unmapped_path,
            },
            threads=1,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "classification": classification_output,
                "missing": missing_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


if __name__ == "__main__":
    run_classification()
