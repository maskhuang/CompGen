"""
Run Liftoff through the Adapter framework.

This script bridges Snakemake with LiftoffAdapter, providing:
- Version checking (AC4)
- Error classification (AC5)
- Audit record generation (AC6)

Source: Story 5.1 Liftoff Adapter

Usage (via Snakemake):
    snakemake.input.reference_gff: Path to reference GFF3 (may be gzipped)
    snakemake.input.reference_fa: Path to reference FASTA (may be gzipped)
    snakemake.input.target_fa: Path to target FASTA (may be gzipped)
    snakemake.output.lifted_gff: Path to lifted GFF3 output
    snakemake.output.unmapped: Path to unmapped features file
    snakemake.output.stats: Path to statistics TSV
    snakemake.output.run_json: Path to .run.json audit file
    snakemake.params.min_coverage: Minimum coverage threshold
    snakemake.params.min_identity: Minimum identity threshold
    snakemake.params.copies: Allow multiple copies flag
    snakemake.params.flank: Flank sequence length
    snakemake.params.out_dir: Output directory
    snakemake.wildcards.reference: Reference species name
    snakemake.wildcards.target: Target species name
    snakemake.threads: Number of threads
    snakemake.config: Full Snakemake config
"""

import csv
import gzip
import shutil
import subprocess
import time
from pathlib import Path
from tempfile import TemporaryDirectory

from workflow.adapters.liftoff import LiftoffAdapter, parse_liftoff_stats
from workflow.adapters.base import AdapterContext
from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.audit import create_and_write_audit


def decompress_if_gzipped(input_path: Path, temp_dir: Path, prefix: str = "file") -> Path:
    """
    Decompress a gzipped file if needed.

    Args:
        input_path: Path to input file (may or may not be gzipped).
        temp_dir: Temporary directory for decompressed file.
        prefix: Prefix for decompressed filename.

    Returns:
        Path to the decompressed file (or original if not gzipped).
    """
    if not str(input_path).endswith('.gz'):
        return input_path

    # Get the base name without .gz to determine file type from actual suffix
    base_name = input_path.name[:-3]  # Remove .gz

    # Determine appropriate extension from the actual file suffix
    if base_name.endswith('.gff3'):
        ext = '.gff3'
    elif base_name.endswith('.gff'):
        ext = '.gff3'  # Normalize to gff3
    elif base_name.endswith(('.fa', '.fasta', '.fna')):
        ext = '.fa'
    else:
        # Fallback: use original extension or .txt
        parts = base_name.rsplit('.', 1)
        ext = f'.{parts[1]}' if len(parts) > 1 else '.txt'

    decompressed_path = temp_dir / f'{prefix}{ext}'
    with gzip.open(input_path, 'rb') as f_in:
        with open(decompressed_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decompressed_path


def write_stats_tsv(stats: dict, output_path: Path, reference: str, target: str) -> None:
    """
    Write Liftoff statistics to TSV file.

    Args:
        stats: Statistics dictionary from parse_liftoff_stats.
        output_path: Path to write TSV file.
        reference: Reference species name.
        target: Target species name.
    """
    fieldnames = [
        "reference", "target", "lifted_genes", "lifted_features",
        "unmapped_genes", "total_genes", "lift_rate"
    ]

    row = {
        "reference": reference,
        "target": target,
        "lifted_genes": stats.get("lifted_genes", 0),
        "lifted_features": stats.get("lifted_features", 0),
        "unmapped_genes": stats.get("unmapped_genes", 0),
        "total_genes": stats.get("total_genes", 0),
        "lift_rate": stats.get("lift_rate", 0.0),
    }

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerow(row)


def run_liftoff_with_adapter():
    """
    Execute Liftoff using the adapter pattern.

    This function:
    1. Validates Liftoff version (AC4)
    2. Decompresses input files if needed
    3. Builds and executes the Liftoff command
    4. Classifies errors if execution fails (AC5)
    5. Generates statistics TSV
    6. Generates .run.json audit record (AC6)
    """
    start_time = time.time()

    # Get Snakemake variables
    ref_gff_path = Path(snakemake.input.reference_gff).resolve()
    ref_fa_path = Path(snakemake.input.reference_fa).resolve()
    target_fa_path = Path(snakemake.input.target_fa).resolve()

    lifted_gff_output = Path(snakemake.output.lifted_gff).resolve()
    unmapped_output = Path(snakemake.output.unmapped).resolve()
    stats_output = Path(snakemake.output.stats).resolve()
    run_json_output = Path(snakemake.output.run_json).resolve()

    reference = snakemake.wildcards.reference
    target = snakemake.wildcards.target
    threads = snakemake.threads

    out_dir = Path(snakemake.params.out_dir)

    # Build config dict for adapter
    config = {
        "liftoff": {
            "min_coverage": snakemake.params.min_coverage,
            "min_identity": snakemake.params.min_identity,
            "copies": snakemake.params.copies,
            "flank": snakemake.params.flank,
        }
    }

    # Initialize adapter
    adapter = LiftoffAdapter()

    # Track execution for audit
    cmd = []
    version = "unknown"
    exit_code = 0
    error_code_str = None
    error_message = None

    try:
        # Step 1: Check version (AC4)
        version = adapter.check_version()
        print(f"Liftoff version: {version}")

        # Step 2: Prepare working directory (use absolute paths)
        out_dir = out_dir.resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        intermediate_dir = (out_dir / "intermediate").resolve()
        intermediate_dir.mkdir(parents=True, exist_ok=True)

        with TemporaryDirectory(dir=out_dir) as temp_dir:
            temp_path = Path(temp_dir)

            # Decompress input files if gzipped
            actual_ref_gff = decompress_if_gzipped(ref_gff_path, temp_path, "reference")
            actual_ref_fa = decompress_if_gzipped(ref_fa_path, temp_path, "reference_genome")
            actual_target_fa = decompress_if_gzipped(target_fa_path, temp_path, "target_genome")

            print(f"Reference GFF: {actual_ref_gff}")
            print(f"Reference FASTA: {actual_ref_fa}")
            print(f"Target FASTA: {actual_target_fa}")

            # Step 3: Build context for adapter
            ctx = AdapterContext(
                inputs={
                    "reference_gff": actual_ref_gff,
                    "reference_fa": actual_ref_fa,
                    "target_fa": actual_target_fa,
                },
                outputs={
                    "lifted_gff": lifted_gff_output,
                    "unmapped": unmapped_output,
                    "intermediate_dir": intermediate_dir,
                },
                config=config,
                wildcards={"reference": reference, "target": target},
                threads=threads,
            )

            # Step 4: Validate inputs
            adapter.validate_inputs(ctx)

            # Step 5: Build and execute command
            cmd = adapter.build_command(ctx)
            print(f"Executing: {' '.join(cmd)}")

            # Set up environment with additional paths for minimap2
            import os
            env = os.environ.copy()
            extra_paths = [
                str(Path.home() / ".local" / "bin"),
                str(Path.home() / "software" / "bin"),
            ]
            env["PATH"] = ":".join(extra_paths) + ":" + env.get("PATH", "")

            result = subprocess.run(
                cmd,
                cwd=out_dir,
                capture_output=True,
                text=True,
                timeout=adapter.timeout_seconds(ctx),
                env=env,
            )

            exit_code = result.returncode
            print(f"Liftoff stdout:\n{result.stdout}")
            if result.stderr:
                print(f"Liftoff stderr:\n{result.stderr}")

            if exit_code != 0:
                # Step 5b: Classify error (AC5)
                error_code, is_retryable = adapter.classify_error(
                    ctx, exit_code, result.stderr
                )
                error_code_str = error_code.value
                error_message = f"Liftoff failed with exit code {exit_code}"

                raise CompGeneError(
                    error_code,
                    error_message,
                    details=result.stderr[:500] if result.stderr else None,
                )

            # Step 6: Verify output files exist
            if not lifted_gff_output.exists():
                raise CompGeneError(
                    ErrorCode.E_OUTPUT_MISSING,
                    f"Liftoff output not found: {lifted_gff_output}",
                )

            # Create unmapped file if it doesn't exist (Liftoff may not create it)
            if not unmapped_output.exists():
                unmapped_output.write_text("")

            # Step 7: Generate statistics
            stats = parse_liftoff_stats(lifted_gff_output, unmapped_output)
            write_stats_tsv(stats, stats_output, reference, target)

            print(f"Liftoff completed: {stats['lifted_genes']} genes lifted "
                  f"({stats['lifted_features']} total features), "
                  f"{stats['unmapped_genes']} genes unmapped, "
                  f"lift rate: {stats['lift_rate']:.2%}")

    except CompGeneError as e:
        error_code_str = e.error_code.value
        error_message = str(e)
        exit_code = e.to_exit_code()
        raise

    except subprocess.TimeoutExpired as e:
        error_code_str = ErrorCode.E_TIMEOUT.value
        timeout_val = e.timeout if e.timeout else "unknown"
        error_message = f"Liftoff timed out after {timeout_val} seconds"
        exit_code = 124  # Standard timeout exit code
        raise CompGeneError(ErrorCode.E_TIMEOUT, error_message)

    finally:
        # Step 8: Write audit record (AC6)
        runtime = time.time() - start_time

        _, audit_path = create_and_write_audit(
            rule="liftoff_map",
            wildcards={"reference": reference, "target": target},
            cmd=cmd,
            tool_version=version,
            input_paths={
                "reference_gff": ref_gff_path,
                "reference_fa": ref_fa_path,
                "target_fa": target_fa_path,
            },
            threads=threads,
            runtime_seconds=runtime,
            exit_code=exit_code,
            meta_dir=run_json_output.parent,
            output_paths={
                "lifted_gff": lifted_gff_output,
                "unmapped": unmapped_output,
                "stats": stats_output,
            } if exit_code == 0 else None,
            error_code=error_code_str,
            error_message=error_message,
        )

        # Move audit file to expected location if different
        if audit_path != run_json_output:
            shutil.move(audit_path, run_json_output)

        print(f"Audit record written to: {run_json_output}")


if __name__ == "__main__":
    run_liftoff_with_adapter()
