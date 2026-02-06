"""
CompGene Orthogroup Utilities Module.

This module provides functionality for parsing OrthoFinder output
and generating structured orthogroup tables for downstream analysis.

Key features:
- Parse OrthoFinder Orthogroups.tsv to long format (one gene per row)
- Calculate per-orthogroup statistics (gene count, species count, single-copy)
- Generate species overlap matrix (shared orthogroup counts)
- Atomic file writing for checkpoint safety

Source: Story 3.2 - Orthogroups 表生成
"""

import csv
import logging
from collections import defaultdict
from pathlib import Path

from workflow.lib.errors import ErrorCode, CompGeneError

logger = logging.getLogger(__name__)

# Module version for audit records
__version__ = "1.0.0"


# =============================================================================
# Parsing
# =============================================================================

def parse_orthogroups_tsv(orthogroups_path: Path) -> list[dict]:
    """
    Parse OrthoFinder Orthogroups.tsv into long-format records.

    OrthoFinder format:
        Orthogroup  species1        species2        species3
        OG0000000   gene1, gene2    gene3           gene4, gene5

    Genes within a cell are comma-separated (`, `). Empty cells mean the
    species has no gene in that orthogroup.

    Args:
        orthogroups_path: Path to OrthoFinder Orthogroups.tsv.

    Returns:
        List of dicts with keys: orthogroup_id, gene_id, species_id.

    Raises:
        CompGeneError: If file does not exist.
    """
    if not orthogroups_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Orthogroups file not found: {orthogroups_path}",
        )

    records: list[dict] = []

    with open(orthogroups_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")

        # Parse header to get species names
        header = next(reader)
        # Header: Orthogroup, species1, species2, ...
        species_names = [col.strip() for col in header[1:]]

        for row in reader:
            if not row or not row[0].strip():
                continue

            orthogroup_id = row[0].strip()
            cells = row[1:]

            for i, cell in enumerate(cells):
                if i >= len(species_names):
                    break
                species_id = species_names[i]
                cell = cell.strip()
                if not cell:
                    continue

                # Split comma-separated gene IDs
                gene_ids = [g.strip() for g in cell.split(",")]
                for gene_id in gene_ids:
                    if gene_id:
                        records.append({
                            "orthogroup_id": orthogroup_id,
                            "gene_id": gene_id,
                            "species_id": species_id,
                        })

    logger.info(
        "Parsed %d gene records from %d orthogroups in %s",
        len(records),
        len({r["orthogroup_id"] for r in records}),
        orthogroups_path.name,
    )
    return records


# =============================================================================
# Statistics
# =============================================================================

def calculate_orthogroup_stats(records: list[dict]) -> list[dict]:
    """
    Calculate per-orthogroup statistics from long-format records.

    For each orthogroup, computes:
    - gene_count: Total number of genes
    - species_count: Number of species present
    - is_single_copy: True if every species has exactly 1 gene AND >=2 species
    - species_list: Comma-separated sorted species names

    Args:
        records: Long-format records from parse_orthogroups_tsv.

    Returns:
        List of stat dicts sorted by orthogroup_id.
    """
    if not records:
        return []

    # Group by orthogroup
    og_genes: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        og_genes[rec["orthogroup_id"]].append(rec)

    stats = []
    for og_id in sorted(og_genes.keys()):
        genes = og_genes[og_id]
        species_set = set()
        species_gene_counts: dict[str, int] = defaultdict(int)

        for g in genes:
            species_set.add(g["species_id"])
            species_gene_counts[g["species_id"]] += 1

        gene_count = len(genes)
        species_count = len(species_set)

        # Single-copy: every species has exactly 1 gene AND at least 2 species
        is_single_copy = (
            species_count >= 2
            and all(c == 1 for c in species_gene_counts.values())
        )

        stats.append({
            "orthogroup_id": og_id,
            "gene_count": gene_count,
            "species_count": species_count,
            "is_single_copy": is_single_copy,
            "species_list": ",".join(sorted(species_set)),
        })

    return stats


def calculate_species_overlap(records: list[dict]) -> dict[str, dict[str, int]]:
    """
    Calculate species overlap matrix: shared orthogroup counts between species pairs.

    The matrix is symmetric. Diagonal values represent the total number of
    orthogroups each species participates in.

    Args:
        records: Long-format records from parse_orthogroups_tsv.

    Returns:
        Nested dict: matrix[species_a][species_b] = shared orthogroup count.
    """
    if not records:
        return {}

    # Build orthogroup → set of species
    og_species: dict[str, set[str]] = defaultdict(set)
    for rec in records:
        og_species[rec["orthogroup_id"]].add(rec["species_id"])

    all_species = sorted({rec["species_id"] for rec in records})

    # Initialize matrix
    matrix: dict[str, dict[str, int]] = {
        sp: {sp2: 0 for sp2 in all_species}
        for sp in all_species
    }

    # Count shared orthogroups
    for species_set in og_species.values():
        species_list = sorted(species_set)
        for i, s1 in enumerate(species_list):
            for s2 in species_list[i:]:
                matrix[s1][s2] += 1
                if s1 != s2:
                    matrix[s2][s1] += 1

    return matrix


# =============================================================================
# File Output
# =============================================================================

def write_orthogroups_long_format(records: list[dict], output_path: Path) -> None:
    """
    Write orthogroups in long format (one gene per row) to TSV.

    Columns: orthogroup_id, gene_id, species_id

    Args:
        records: Long-format records from parse_orthogroups_tsv.
        output_path: Destination file path.
    """
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(tmp_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["orthogroup_id", "gene_id", "species_id"])
            for rec in records:
                writer.writerow([rec["orthogroup_id"], rec["gene_id"], rec["species_id"]])
        tmp_path.rename(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    logger.info("Wrote %d records to %s", len(records), output_path)


def write_orthogroup_stats(stats: list[dict], output_path: Path) -> None:
    """
    Write orthogroup statistics to TSV.

    Columns: orthogroup_id, gene_count, species_count, is_single_copy, species_list

    Args:
        stats: Stats from calculate_orthogroup_stats.
        output_path: Destination file path.
    """
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(tmp_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["orthogroup_id", "gene_count", "species_count", "is_single_copy", "species_list"])
            for s in stats:
                writer.writerow([
                    s["orthogroup_id"],
                    s["gene_count"],
                    s["species_count"],
                    str(s["is_single_copy"]).lower(),
                    s["species_list"],
                ])
        tmp_path.rename(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    logger.info("Wrote %d orthogroup stats to %s", len(stats), output_path)


def write_species_overlap(matrix: dict[str, dict[str, int]], output_path: Path) -> None:
    """
    Write species overlap matrix to TSV.

    First column is 'species', followed by one column per species.
    Values are shared orthogroup counts.

    Args:
        matrix: Overlap matrix from calculate_species_overlap.
        output_path: Destination file path.
    """
    if not matrix:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("species\n")
        return

    species = sorted(matrix.keys())
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(tmp_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["species"] + species)
            for sp in species:
                row = [sp] + [str(matrix[sp][sp2]) for sp2 in species]
                writer.writerow(row)
        tmp_path.rename(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    logger.info("Wrote %dx%d species overlap matrix to %s", len(species), len(species), output_path)
