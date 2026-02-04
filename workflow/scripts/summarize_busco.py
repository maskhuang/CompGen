"""
Summarize BUSCO results across all species.

This script aggregates BUSCO short_summary.txt files from multiple species
into a single TSV file for cross-species comparison.

Source: Story 2.4 BUSCO QC (AC3)

Usage (via Snakemake):
    snakemake.input.summaries: List of short_summary.txt file paths
    snakemake.output.summary: Path to output TSV file

Output columns:
    - species: Species identifier
    - complete_pct: Percentage of complete BUSCOs
    - single_copy_pct: Percentage of single-copy BUSCOs
    - duplicated_pct: Percentage of duplicated BUSCOs
    - fragmented_pct: Percentage of fragmented BUSCOs
    - missing_pct: Percentage of missing BUSCOs
    - total: Total BUSCO groups searched
    - lineage: Lineage dataset used
    - busco_version: BUSCO version used
"""

import csv
from pathlib import Path

# Import shared BUSCO parsing patterns from adapter to avoid duplication
from workflow.adapters.busco import (
    BUSCO_SUMMARY_PATTERN,
    LINEAGE_PATTERN,
    VERSION_PATTERN,
)


# =============================================================================
# Parsing functions
# =============================================================================

def parse_summary(path: Path) -> dict:
    """
    Parse a single BUSCO short_summary.txt file.

    Args:
        path: Path to short_summary.txt file.

    Returns:
        Dictionary with parsed statistics.

    Raises:
        ValueError: If BUSCO notation cannot be parsed.
    """
    content = path.read_text()

    # Parse statistics from BUSCO notation
    match = BUSCO_SUMMARY_PATTERN.search(content)
    if not match:
        raise ValueError(f"Cannot parse BUSCO summary: {path}")

    # Parse lineage
    lineage_match = LINEAGE_PATTERN.search(content)
    lineage = lineage_match.group(1) if lineage_match else "unknown"

    # Parse version
    version_match = VERSION_PATTERN.search(content)
    version = version_match.group(1) if version_match else "unknown"

    # Extract species name from path
    # Expected path: results/qc/{species}/busco/short_summary.txt
    species = path.parent.parent.name

    return {
        "species": species,
        "complete_pct": float(match.group(1)),
        "single_copy_pct": float(match.group(2)),
        "duplicated_pct": float(match.group(3)),
        "fragmented_pct": float(match.group(4)),
        "missing_pct": float(match.group(5)),
        "total": int(match.group(6)),
        "lineage": lineage,
        "busco_version": version,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    """Main entry point for Snakemake script."""
    # Get inputs and outputs from Snakemake
    summaries = snakemake.input.summaries
    output_path = Path(snakemake.output.summary)

    # Parse all summary files
    results = []
    for summary_path in summaries:
        try:
            result = parse_summary(Path(summary_path))
            results.append(result)
        except Exception as e:
            print(f"Warning: Failed to parse {summary_path}: {e}")

    # Sort by species name
    results.sort(key=lambda x: x["species"])

    # Atomic write: write to .tmp then rename
    temp_path = output_path.with_suffix('.tmp')

    fieldnames = [
        "species",
        "complete_pct",
        "single_copy_pct",
        "duplicated_pct",
        "fragmented_pct",
        "missing_pct",
        "total",
        "lineage",
        "busco_version",
    ]

    with open(temp_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    # Atomic rename
    temp_path.rename(output_path)

    print(f"BUSCO summary written to {output_path} ({len(results)} species)")


if __name__ == "__main__":
    main()
