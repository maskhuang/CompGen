"""
Build structured orthogroup tables from OrthoFinder output.

This script bridges Snakemake with lib/orthogroup_utils.py, generating:
- orthogroups.tsv: Long-format gene-to-orthogroup mapping
- orthogroup_stats.tsv: Per-orthogroup statistics
- species_overlap.tsv: Pairwise species overlap matrix
- .run.json: Audit record

Source: Story 3.2 - Orthogroups 表生成

Usage (via Snakemake):
    snakemake.input.orthogroups_raw: OrthoFinder Orthogroups.tsv
    snakemake.output.orthogroups: Long-format orthogroups.tsv
    snakemake.output.stats: orthogroup_stats.tsv
    snakemake.output.overlap: species_overlap.tsv
    snakemake.output.run_json: .run.json audit file
"""

import shutil
import time
from pathlib import Path

from workflow.lib.orthogroup_utils import (
    __version__,
    parse_orthogroups_tsv,
    calculate_orthogroup_stats,
    calculate_species_overlap,
    write_orthogroups_long_format,
    write_orthogroup_stats,
    write_species_overlap,
)
from workflow.lib.audit import create_and_write_audit
from workflow.lib.errors import CompGeneError


def main():
    start_time = time.time()

    # Get Snakemake variables
    orthogroups_raw = Path(snakemake.input.orthogroups_raw)
    orthogroups_output = Path(snakemake.output.orthogroups)
    stats_output = Path(snakemake.output.stats)
    overlap_output = Path(snakemake.output.overlap)
    run_json_output = Path(snakemake.output.run_json)

    # Track execution for audit
    exit_code = 0
    error_code_str = None
    error_message = None

    try:
        # 1. Parse OrthoFinder output
        print(f"Parsing OrthoFinder orthogroups from {orthogroups_raw}")
        records = parse_orthogroups_tsv(orthogroups_raw)
        print(f"Parsed {len(records)} gene records")

        # 2. Calculate statistics
        stats = calculate_orthogroup_stats(records)
        overlap = calculate_species_overlap(records)
        num_species = len(overlap) if overlap else 0
        single_copy = sum(1 for s in stats if s["is_single_copy"])
        print(f"Statistics: {len(stats)} orthogroups, {num_species} species, "
              f"{single_copy} single-copy orthogroups")

        # 3. Write output files
        write_orthogroups_long_format(records, orthogroups_output)
        print(f"Long-format table written to: {orthogroups_output}")

        write_orthogroup_stats(stats, stats_output)
        print(f"Orthogroup stats written to: {stats_output}")

        write_species_overlap(overlap, overlap_output)
        print(f"Species overlap matrix written to: {overlap_output}")

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
        # 4. Generate audit record
        runtime = time.time() - start_time

        run_json_output.parent.mkdir(parents=True, exist_ok=True)

        cmd = ["python", "build_orthogroup_tables.py"]

        _, audit_path = create_and_write_audit(
            rule="orthology_parse_orthogroups",
            wildcards={},
            cmd=cmd,
            tool_version=__version__,
            input_paths={"orthogroups_raw": orthogroups_raw},
            threads=1,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "orthogroups": orthogroups_output,
                "stats": stats_output,
                "overlap": overlap_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


if __name__ == "__main__":
    main()
