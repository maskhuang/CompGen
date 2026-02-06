"""
Build presence/absence and copy number matrices from orthogroups.

This script bridges Snakemake with lib/presence_absence.py and
lib/orthogroup_utils.py, generating:
- presence_absence.tsv: Binary presence/absence matrix
- copy_number.tsv: Gene copy number matrix
- .run.json: Audit record

Source: Story 6A.1 - Presence/Absence 矩阵生成

Usage (via Snakemake):
    snakemake.input.orthogroups: Long-format orthogroups.tsv
    snakemake.output.presence_absence: Binary matrix
    snakemake.output.copy_number: Count matrix
    snakemake.output.run_json: .run.json audit file
"""

import shutil
import time
from pathlib import Path

from workflow.lib.presence_absence import (
    __version__,
    build_copy_number_matrix,
    get_sorted_species,
    read_orthogroups_long_format,
    write_matrix_tsv,
)
from workflow.lib.audit import create_and_write_audit
from workflow.lib.errors import CompGeneError


def main():
    start_time = time.time()

    # Get Snakemake variables
    orthogroups_input = Path(snakemake.input.orthogroups)
    pa_output = Path(snakemake.output.presence_absence)
    cn_output = Path(snakemake.output.copy_number)
    run_json_output = Path(snakemake.output.run_json)

    # Track execution for audit
    exit_code = 0
    error_code_str = None
    error_message = None

    try:
        # 1. Parse orthogroups long-format table (Story 3.2 output)
        print(f"Reading orthogroups from {orthogroups_input}")
        records = read_orthogroups_long_format(orthogroups_input)
        species = get_sorted_species(records)
        print(f"Parsed {len(records)} gene records, {len(species)} species")

        # 2. Build matrices (CN first, PA derived from CN to avoid double computation)
        cn_matrix = build_copy_number_matrix(records)
        pa_matrix = {
            og: {sp: (1 if v > 0 else 0) for sp, v in row.items()}
            for og, row in cn_matrix.items()
        }
        print(f"Generated matrices: {len(cn_matrix)} orthogroups x {len(species)} species")

        # 3. Write output files
        write_matrix_tsv(pa_matrix, species, pa_output)
        print(f"Presence/absence matrix written to: {pa_output}")

        write_matrix_tsv(cn_matrix, species, cn_output)
        print(f"Copy number matrix written to: {cn_output}")

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

        _, audit_path = create_and_write_audit(
            rule="matrices_generate",
            wildcards={},
            cmd=["python", "build_matrices.py"],
            tool_version=__version__,
            input_paths={"orthogroups": orthogroups_input},
            threads=1,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "presence_absence": pa_output,
                "copy_number": cn_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


main()
