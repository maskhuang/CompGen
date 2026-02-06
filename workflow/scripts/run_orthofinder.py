"""
Run OrthoFinder through the Adapter framework.

This script bridges Snakemake with OrthoFinderAdapter, providing:
- Version checking (AC4)
- Error classification (AC5)
- Audit record generation (AC6)
- Output directory management (OrthoFinder creates Results_compgene/ subdirectory)

Source: Story 3.1 OrthoFinder Adapter

Usage (via Snakemake):
    snakemake.input.proteins_dir: Path to directory with .fa files per species
    snakemake.input.done: Sentinel file from prepare step
    snakemake.output.orthogroups: Path to Orthogroups.tsv
    snakemake.output.gene_count: Path to Orthogroups.GeneCount.tsv
    snakemake.output.species_tree: Path to SpeciesTree_rooted.txt
    snakemake.output.run_json: Path to .run.json audit file
    snakemake.params.results_dir: Base results directory for orthology
    snakemake.threads: Number of threads
    snakemake.config: Full Snakemake config
"""

import os
import shutil
import subprocess
import time
from pathlib import Path

from workflow.adapters.orthofinder import OrthoFinderAdapter
from workflow.adapters.base import AdapterContext
from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.audit import create_and_write_audit


def find_orthofinder_results(output_dir: Path) -> Path:
    """
    Locate the OrthoFinder results directory.

    OrthoFinder writes results to:
    - With -o flag: {output_dir}/Results_compgene/
    - Without -o: {input_dir}/OrthoFinder/Results_compgene/

    Args:
        output_dir: The -o output directory passed to OrthoFinder.

    Returns:
        Path to the Results_compgene directory.

    Raises:
        CompGeneError: If results directory cannot be found.
    """
    # With -n compgene, results are in Results_compgene/
    results_dir = output_dir / "Results_compgene"
    if results_dir.is_dir():
        return results_dir

    # Fallback: search for any Results_* directory
    for d in sorted(output_dir.iterdir()):
        if d.is_dir() and d.name.startswith("Results_"):
            return d

    raise CompGeneError(
        ErrorCode.E_OUTPUT_MISSING,
        f"OrthoFinder results directory not found in {output_dir}",
    )


def link_results(src_dir: Path, dest_dir: Path) -> None:
    """
    Move OrthoFinder output subdirectories to the expected Snakemake output paths.

    OrthoFinder creates its own directory structure under Results_compgene/.
    This function moves the relevant subdirectories to the expected locations.

    Args:
        src_dir: OrthoFinder Results_compgene/ directory.
        dest_dir: Expected output directory (results/orthology/).
    """
    subdirs_to_move = [
        "Orthogroups",
        "Species_Tree",
        "Gene_Trees",
        "Resolved_Gene_Trees",
        "Phylogenetic_Hierarchical_Orthogroups",
        "Comparative_Genomics_Statistics",
    ]

    for subdir in subdirs_to_move:
        src = src_dir / subdir
        dest = dest_dir / subdir
        if src.is_dir():
            if dest.exists():
                shutil.rmtree(str(dest))
            shutil.move(str(src), str(dest))


def main():
    adapter = OrthoFinderAdapter()

    # 1. Version check (AC4)
    version = adapter.check_version()

    # 2. Build context from snakemake variables
    proteins_dir = Path(snakemake.input.proteins_dir).resolve()
    results_dir = Path(snakemake.params.results_dir).resolve()

    # OrthoFinder -o writes to a new directory; use a staging area
    of_output_dir = results_dir / "_orthofinder_run"
    of_output_dir.mkdir(parents=True, exist_ok=True)

    ctx = AdapterContext(
        inputs={"proteins_dir": proteins_dir},
        outputs={
            "results_dir": results_dir,
            "output_dir": of_output_dir,
        },
        config=snakemake.config,
        wildcards=dict(snakemake.wildcards) if snakemake.wildcards else {},
        threads=snakemake.threads,
    )

    # 3. Validate inputs (AC5)
    adapter.validate_inputs(ctx)

    # 4. Build and execute command
    cmd = adapter.build_command(ctx)
    start_time = time.time()

    log_path = Path(snakemake.log[0]).resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        timeout = adapter.timeout_seconds(ctx)

        # Redirect stdout/stderr to log file to avoid pipe buffer deadlock
        # (lesson from Story 5.1 Liftoff Adapter)
        with open(str(log_path), 'w') as log_fh:
            proc = subprocess.run(
                cmd,
                stdout=log_fh,
                stderr=log_fh,
                timeout=timeout,
            )

        runtime = time.time() - start_time

        if proc.returncode != 0:
            log_content = log_path.read_text() if log_path.exists() else ""
            error_code, retryable = adapter.classify_error(ctx, proc.returncode, log_content)
            raise CompGeneError(
                error_code,
                f"OrthoFinder failed (exit {proc.returncode}). Check log: {log_path}",
            )

    except subprocess.TimeoutExpired as e:
        runtime = time.time() - start_time
        raise CompGeneError(
            ErrorCode.E_TIMEOUT,
            f"OrthoFinder timed out after {e.timeout}s. "
            f"Increase orthofinder.timeout in config or reduce species count.",
        )

    # 5. Move OrthoFinder output to expected paths
    of_results = find_orthofinder_results(of_output_dir)
    link_results(of_results, results_dir)

    # Clean up staging directory
    shutil.rmtree(str(of_output_dir), ignore_errors=True)

    # 6. Update context to point to final locations for parse_outputs
    ctx_final = AdapterContext(
        inputs={"proteins_dir": proteins_dir},
        outputs={"results_dir": results_dir},
        config=snakemake.config,
        wildcards=dict(snakemake.wildcards) if snakemake.wildcards else {},
        threads=snakemake.threads,
    )

    # 7. Parse outputs
    result = adapter.parse_outputs(ctx_final)

    # 8. Generate audit record (AC6)
    audit_path = Path(snakemake.output.run_json).resolve()
    audit_path.parent.mkdir(parents=True, exist_ok=True)

    create_and_write_audit(
        rule="orthology_infer",
        wildcards=dict(snakemake.wildcards) if snakemake.wildcards else {},
        cmd=cmd,
        tool_version=version,
        input_paths={"proteins_dir": proteins_dir},
        threads=snakemake.threads,
        runtime_seconds=runtime,
        exit_code=0,
        meta_dir=audit_path.parent.parent,
        output_paths={
            "orthogroups": Path(snakemake.output.orthogroups),
            "gene_count": Path(snakemake.output.gene_count),
            "species_tree": Path(snakemake.output.species_tree),
        },
    )


if __name__ == "__main__":
    main()
