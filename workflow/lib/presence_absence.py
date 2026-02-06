"""
CompGene Presence/Absence Matrix Module.

This module provides functionality for generating presence/absence
and copy number matrices from orthogroup records.

Key features:
- Read long-format orthogroup TSV (orthogroup_id, gene_id, species_id)
- Build binary presence/absence matrix from long-format orthogroup records
- Build copy number matrix showing gene counts per species per orthogroup
- Atomic file writing for checkpoint safety
- Sorted species columns and orthogroup rows for reproducibility

Source: Story 6A.1 - Presence/Absence 矩阵生成
"""

import csv
import logging
from collections import defaultdict
from pathlib import Path

from workflow.lib.errors import CompGeneError, ErrorCode

logger = logging.getLogger(__name__)

# Module version for audit records
__version__ = "1.0.0"


# =============================================================================
# Input Reading
# =============================================================================

def read_orthogroups_long_format(input_path: Path) -> list[dict]:
    """
    Read long-format orthogroups TSV into records.

    Expected format (output of Story 3.2 orthology_parse_orthogroups):
        orthogroup_id  gene_id    species_id
        OG0000000      gene1      species1
        OG0000000      gene2      species1

    This is NOT the OrthoFinder raw format. Use parse_orthogroups_tsv()
    from orthogroup_utils.py for the raw OrthoFinder Orthogroups.tsv.

    Args:
        input_path: Path to long-format orthogroups.tsv.

    Returns:
        List of dicts with keys: orthogroup_id, gene_id, species_id.

    Raises:
        CompGeneError: E_INPUT_MISSING if file does not exist.
        CompGeneError: E_INPUT_FORMAT if header columns are missing or unexpected.
    """
    if not input_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Orthogroups file not found: {input_path}",
        )

    records: list[dict] = []

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        # Validate header
        expected_cols = {"orthogroup_id", "gene_id", "species_id"}
        if reader.fieldnames is None or not expected_cols.issubset(set(reader.fieldnames)):
            raise CompGeneError(
                ErrorCode.E_INPUT_FORMAT,
                f"Invalid header in {input_path}. "
                f"Expected columns: {expected_cols}, got: {reader.fieldnames}",
            )

        for row in reader:
            og_id = row["orthogroup_id"].strip()
            gene_id = row["gene_id"].strip()
            species_id = row["species_id"].strip()
            if og_id and gene_id and species_id:
                records.append({
                    "orthogroup_id": og_id,
                    "gene_id": gene_id,
                    "species_id": species_id,
                })

    logger.info(
        "Read %d gene records from %d orthogroups in %s",
        len(records),
        len({r["orthogroup_id"] for r in records}),
        input_path.name,
    )
    return records


# =============================================================================
# Species Extraction
# =============================================================================

def get_sorted_species(records: list[dict]) -> list[str]:
    """
    Extract and sort unique species names from long-format records.

    Args:
        records: Long-format records with 'species_id' key.

    Returns:
        Alphabetically sorted list of unique species names.
    """
    if not records:
        return []
    return sorted({rec["species_id"] for rec in records})


# =============================================================================
# Matrix Generation
# =============================================================================

def build_copy_number_matrix(records: list[dict]) -> dict[str, dict[str, int]]:
    """
    Build copy number matrix from long-format orthogroup records.

    Each cell contains the number of genes a species has in that orthogroup.
    Species absent from an orthogroup get value 0.

    Args:
        records: Long-format records with orthogroup_id, gene_id, species_id.

    Returns:
        Nested dict: matrix[orthogroup_id][species_id] = gene_count.
        All species appear as columns in every orthogroup row.
    """
    if not records:
        return {}

    all_species = get_sorted_species(records)

    # Count genes per (orthogroup, species) pair
    counts: dict[str, dict[str, int]] = defaultdict(lambda: {sp: 0 for sp in all_species})
    for rec in records:
        counts[rec["orthogroup_id"]][rec["species_id"]] += 1

    return dict(counts)


def build_presence_absence_matrix(records: list[dict]) -> dict[str, dict[str, int]]:
    """
    Build binary presence/absence matrix from long-format orthogroup records.

    Each cell is 1 if the species has at least one gene in the orthogroup,
    0 otherwise. Multi-copy genes still yield 1.

    Args:
        records: Long-format records with orthogroup_id, gene_id, species_id.

    Returns:
        Nested dict: matrix[orthogroup_id][species_id] = 0 or 1.
        All species appear as columns in every orthogroup row.
    """
    cn_matrix = build_copy_number_matrix(records)

    pa_matrix: dict[str, dict[str, int]] = {}
    for og_id, row in cn_matrix.items():
        pa_matrix[og_id] = {sp: (1 if count > 0 else 0) for sp, count in row.items()}

    return pa_matrix


# =============================================================================
# File Output
# =============================================================================

def write_matrix_tsv(
    matrix: dict[str, dict[str, int]],
    species_list: list[str],
    output_path: Path,
) -> None:
    """
    Write a matrix (presence/absence or copy number) to TSV.

    Format:
        orthogroup_id  species1  species2  species3
        OG0000000      1         1         0

    Rows are sorted by orthogroup_id. Columns follow species_list order.
    Uses atomic write (write to .tmp, then rename) for checkpoint safety.

    Args:
        matrix: Nested dict matrix[orthogroup_id][species_id] = value.
        species_list: Ordered list of species for column headers.
        output_path: Destination file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    try:
        with open(tmp_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["orthogroup_id"] + species_list)
            for og_id in sorted(matrix.keys()):
                row = [og_id] + [str(matrix[og_id].get(sp, 0)) for sp in species_list]
                writer.writerow(row)
        tmp_path.rename(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    logger.info(
        "Wrote %dx%d matrix to %s",
        len(matrix),
        len(species_list),
        output_path,
    )
