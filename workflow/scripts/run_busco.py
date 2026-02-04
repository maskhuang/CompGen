"""
Run BUSCO through the Adapter framework.

This script bridges Snakemake with BuscoAdapter, providing:
- Version checking (AC4)
- Error classification (AC5)
- Audit record generation (AC6)

Source: Story 2.4 BUSCO QC

Usage (via Snakemake):
    snakemake.input.proteins: Path to protein sequences (gzipped or plain)
    snakemake.output.summary: Path to short_summary.txt
    snakemake.output.full_table: Path to full_table.tsv
    snakemake.output.missing_list: Path to missing_busco_list.tsv
    snakemake.output.run_json: Path to .run.json audit file
    snakemake.params.lineage: BUSCO lineage database
    snakemake.params.mode: BUSCO mode (proteins/genome)
    snakemake.params.download_path: Lineage database path (optional)
    snakemake.params.offline: Offline mode flag
    snakemake.params.out_dir: BUSCO output directory
    snakemake.wildcards.species: Species identifier
    snakemake.threads: Number of threads
    snakemake.config: Full Snakemake config
"""

import gzip
import shutil
import subprocess
import time
from pathlib import Path
from tempfile import TemporaryDirectory

from workflow.adapters.busco import BuscoAdapter
from workflow.adapters.base import AdapterContext, RunResult
from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.audit import create_and_write_audit


def decompress_if_gzipped(input_path: Path, temp_dir: Path) -> Path:
    """
    Decompress a gzipped file if needed.

    Args:
        input_path: Path to input file (may or may not be gzipped).
        temp_dir: Temporary directory for decompressed file.

    Returns:
        Path to the decompressed file (or original if not gzipped).
    """
    if str(input_path).endswith('.gz'):
        decompressed_path = temp_dir / 'proteins.fa'
        with gzip.open(input_path, 'rb') as f_in:
            with open(decompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return decompressed_path
    return input_path


def run_busco_with_adapter():
    """
    Execute BUSCO using the adapter pattern.

    This function:
    1. Validates BUSCO version (AC4)
    2. Builds and executes the BUSCO command
    3. Classifies errors if execution fails (AC5)
    4. Generates .run.json audit record (AC6)
    """
    start_time = time.time()

    # Get Snakemake variables
    proteins_path = Path(snakemake.input.proteins)
    summary_output = Path(snakemake.output.summary)
    full_table_output = Path(snakemake.output.full_table)
    missing_list_output = Path(snakemake.output.missing_list)
    run_json_output = Path(snakemake.output.run_json)

    species = snakemake.wildcards.species
    threads = snakemake.threads

    # Build config dict for adapter
    config = {
        "busco": {
            "lineage": snakemake.params.lineage,
            "mode": snakemake.params.mode,
            "download_path": snakemake.params.download_path or "",
            "offline": snakemake.params.offline,
        }
    }

    out_dir = Path(snakemake.params.out_dir)

    # Initialize adapter
    adapter = BuscoAdapter()

    # Track execution for audit
    cmd = []
    version = "unknown"
    exit_code = 0
    error_code_str = None
    error_message = None

    try:
        # Step 1: Check version (AC4)
        version = adapter.check_version()
        print(f"BUSCO version: {version}")

        # Step 2: Prepare working directory
        out_dir.mkdir(parents=True, exist_ok=True)

        with TemporaryDirectory(dir=out_dir) as temp_dir:
            temp_path = Path(temp_dir)

            # Decompress protein file if gzipped
            actual_proteins = decompress_if_gzipped(proteins_path, temp_path)
            print(f"Using protein file: {actual_proteins}")

            # Step 3: Build context for adapter
            ctx = AdapterContext(
                inputs={"proteins": actual_proteins},
                outputs={
                    "summary": summary_output,
                    "full_table": full_table_output,
                    "missing_list": missing_list_output,
                },
                config=config,
                wildcards={"species": species},
                threads=threads,
            )

            # Step 4: Validate inputs
            adapter.validate_inputs(ctx)

            # Step 5: Build and execute command
            # Note: BUSCO outputs to current directory, so we need to run from out_dir
            cmd = adapter.build_command(ctx)
            print(f"Executing: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                cwd=out_dir,
                capture_output=True,
                text=True,
                timeout=adapter.timeout_seconds(ctx),
            )

            exit_code = result.returncode

            if exit_code != 0:
                # Step 5b: Classify error (AC5)
                error_code, is_retryable = adapter.classify_error(
                    ctx, exit_code, result.stderr
                )
                error_code_str = error_code.value
                error_message = f"BUSCO failed with exit code {exit_code}"

                raise CompGeneError(
                    error_code,
                    error_message,
                    details=result.stderr[:500] if result.stderr else None,
                )

            # Step 6: Locate and copy output files
            # BUSCO creates: run_{species}/short_summary.txt, run_{species}/full_table.tsv, etc.
            busco_run_dir = out_dir / species

            # Try common BUSCO output locations
            possible_run_dirs = [
                out_dir / f"run_{species}",
                out_dir / species,
            ]

            actual_run_dir = None
            for run_dir in possible_run_dirs:
                if run_dir.exists():
                    actual_run_dir = run_dir
                    break

            if actual_run_dir is None:
                raise CompGeneError(
                    ErrorCode.E_OUTPUT_MISSING,
                    f"BUSCO output directory not found. Tried: {possible_run_dirs}",
                )

            # Copy short_summary.txt
            short_summary_src = actual_run_dir / "short_summary.txt"
            if not short_summary_src.exists():
                # Try nested structure: run_{species}/short_summary.txt
                nested = out_dir / f"run_{species}" / "short_summary.txt"
                if nested.exists():
                    short_summary_src = nested

            if short_summary_src.exists():
                shutil.copy2(short_summary_src, summary_output)
            else:
                raise CompGeneError(
                    ErrorCode.E_OUTPUT_MISSING,
                    f"short_summary.txt not found in {actual_run_dir}",
                )

            # Copy full_table.tsv
            full_table_src = actual_run_dir / "full_table.tsv"
            if full_table_src.exists():
                shutil.copy2(full_table_src, full_table_output)
            else:
                raise CompGeneError(
                    ErrorCode.E_OUTPUT_MISSING,
                    f"full_table.tsv not found in {actual_run_dir}",
                )

            # Copy missing_busco_list.tsv
            missing_list_src = actual_run_dir / "missing_busco_list.tsv"
            if missing_list_src.exists():
                shutil.copy2(missing_list_src, missing_list_output)
            else:
                raise CompGeneError(
                    ErrorCode.E_OUTPUT_MISSING,
                    f"missing_busco_list.tsv not found in {actual_run_dir}",
                )

            print(f"BUSCO completed successfully for {species}")

    except CompGeneError as e:
        error_code_str = e.error_code.value
        error_message = str(e)
        exit_code = e.to_exit_code()
        raise

    except subprocess.TimeoutExpired:
        error_code_str = ErrorCode.E_TIMEOUT.value
        error_message = f"BUSCO timed out after {adapter.timeout_seconds(ctx)} seconds"
        exit_code = 124  # Standard timeout exit code
        raise CompGeneError(ErrorCode.E_TIMEOUT, error_message)

    finally:
        # Step 7: Write audit record (AC6)
        runtime = time.time() - start_time

        _, audit_path = create_and_write_audit(
            rule="qc_busco",
            wildcards={"species": species},
            cmd=cmd,
            tool_version=version,
            input_paths={"proteins": proteins_path},
            threads=threads,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "summary": summary_output,
                "full_table": full_table_output,
                "missing_list": missing_list_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


if __name__ == "__main__":
    run_busco_with_adapter()
