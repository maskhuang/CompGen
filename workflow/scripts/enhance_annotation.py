"""
Enhance lifted annotation with Liftoff classification attributes.

This script bridges Snakemake with the annotation enhancement module, providing:
- GFF3 attribute enhancement with liftoff_coverage, liftoff_identity, liftoff_status
- Filtering by classification status (present only, or including uncertain)
- Audit record generation

Source: Story 5.3 - 迁移注释输出

Usage (via Snakemake):
    snakemake.input.lifted_gff: Path to Liftoff lifted_annotation.gff3
    snakemake.input.classification: Path to gene_classification.tsv
    snakemake.output.enhanced_gff: Path to output enhanced GFF3
    snakemake.output.run_json: Path to .run.json audit file
    snakemake.params.include_uncertain: Whether to include uncertain genes
    snakemake.wildcards.reference: Reference species name
    snakemake.wildcards.target: Target species name
"""

import shutil
import time
from pathlib import Path

from workflow.lib.annotation_enhance import (
    __version__ as annotation_enhance_version,
    enhance_gff_with_liftoff_attrs,
)
from workflow.lib.audit import create_and_write_audit
from workflow.lib.errors import CompGeneError


# =============================================================================
# Main Execution
# =============================================================================

def run_enhancement():
    """
    Execute GFF3 enhancement from Liftoff outputs.

    This function:
    1. Reads Liftoff lifted GFF and classification results
    2. Enhances GFF with liftoff_* attributes
    3. Filters to only present (and optionally uncertain) genes
    4. Generates .run.json audit record
    """
    start_time = time.time()

    # Get Snakemake variables
    lifted_gff_path = Path(snakemake.input.lifted_gff)
    classification_path = Path(snakemake.input.classification)

    enhanced_gff_output = Path(snakemake.output.enhanced_gff)
    run_json_output = Path(snakemake.output.run_json)

    reference = snakemake.wildcards.reference
    target = snakemake.wildcards.target

    include_uncertain = snakemake.params.include_uncertain

    # Track execution for audit
    exit_code = 0
    error_code_str = None
    error_message = None
    stats = {}

    try:
        # Step 1: Enhance GFF with liftoff attributes
        print(f"Enhancing GFF from {lifted_gff_path}")
        print(f"Classification source: {classification_path}")
        print(f"Include uncertain: {include_uncertain}")

        stats = enhance_gff_with_liftoff_attrs(
            gff_path=lifted_gff_path,
            classification_path=classification_path,
            output_path=enhanced_gff_output,
            reference_species=reference,
            include_uncertain=include_uncertain,
        )

        # Step 2: Log results
        print(f"Enhancement results: "
              f"{stats['input_genes']} input genes, "
              f"{stats['output_genes']} output genes, "
              f"{stats['filtered_genes']} filtered, "
              f"{stats['unclassified_genes']} unclassified")
        print(f"Enhanced GFF written to: {enhanced_gff_output}")

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
        # Step 3: Write audit record
        runtime = time.time() - start_time

        # Ensure meta dir exists
        run_json_output.parent.mkdir(parents=True, exist_ok=True)

        # Build command with parameters for audit traceability
        cmd = [
            "enhance_annotation.py",
            f"--reference={reference}",
            f"--target={target}",
            f"--include-uncertain={include_uncertain}",
        ]

        _, audit_path = create_and_write_audit(
            rule="liftoff_enhance",
            wildcards={"reference": reference, "target": target},
            cmd=cmd,
            tool_version=annotation_enhance_version,
            input_paths={
                "lifted_gff": lifted_gff_path,
                "classification": classification_path,
            },
            threads=1,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "enhanced_gff": enhanced_gff_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
            extra_metadata={
                "stats": stats if stats else None,
            },
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


if __name__ == "__main__":
    run_enhancement()
