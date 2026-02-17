"""
Build multi-species comparison from presence/absence matrix.

This script bridges Snakemake with lib/presence_absence.py, generating:
- comparison.tsv: Per-orthogroup sharing classification
- comparison_summary.tsv: Aggregate sharing statistics
- .run.json: Audit record

Source: Story 6A.2 - 多物种比较

Usage (via Snakemake):
    snakemake.input.presence_absence: PA matrix TSV
    snakemake.output.comparison: Classification TSV
    snakemake.output.summary: Summary TSV
    snakemake.output.run_json: .run.json audit file
"""

import shutil
import time
from pathlib import Path

from workflow.lib.presence_absence import (
    __version__,
    classify_orthogroup_sharing,
    read_presence_absence_matrix,
    summarize_sharing_counts,
    write_comparison_tsv,
    write_summary_tsv,
)
from workflow.lib.audit import create_and_write_audit
from workflow.lib.errors import CompGeneError


def main():
    start_time = time.time()

    # Get Snakemake variables
    pa_input = Path(snakemake.input.presence_absence)
    comparison_output = Path(snakemake.output.comparison)
    summary_output = Path(snakemake.output.summary)
    run_json_output = Path(snakemake.output.run_json)

    # Track execution for audit
    exit_code = 0
    error_code_str = None
    error_message = None

    try:
        # 1. Read presence/absence matrix
        print(f"Reading PA matrix from {pa_input}")
        matrix, species_list = read_presence_absence_matrix(pa_input)
        print(f"Loaded {len(matrix)} orthogroups x {len(species_list)} species")

        # 2. Classify orthogroup sharing patterns
        classifications = classify_orthogroup_sharing(matrix, species_list)
        print(f"Classified {len(classifications)} orthogroups")

        # 3. Generate summary statistics
        summary_rows = summarize_sharing_counts(classifications, species_list)

        # 4. Write output files
        write_comparison_tsv(classifications, comparison_output)
        print(f"Comparison written to: {comparison_output}")

        write_summary_tsv(summary_rows, summary_output)
        print(f"Summary written to: {summary_output}")

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
        # 5. Generate audit record
        runtime = time.time() - start_time

        run_json_output.parent.mkdir(parents=True, exist_ok=True)

        _, audit_path = create_and_write_audit(
            rule="matrices_compare",
            wildcards={},
            cmd=["python", "build_comparison.py"],
            tool_version=__version__,
            input_paths={"presence_absence": pa_input},
            threads=1,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "comparison": comparison_output,
                "summary": summary_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


main()
