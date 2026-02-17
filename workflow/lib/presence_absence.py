"""
CompGene Presence/Absence Matrix Module.

This module provides functionality for generating presence/absence
and copy number matrices from orthogroup records, and for comparing
multi-species sharing patterns.

Key features:
- Read long-format orthogroup TSV (orthogroup_id, gene_id, species_id)
- Build binary presence/absence matrix from long-format orthogroup records
- Build copy number matrix showing gene counts per species per orthogroup
- Read PA matrix and classify orthogroups by sharing pattern
- Summarize sharing statistics across species
- Atomic file writing for checkpoint safety
- Sorted species columns and orthogroup rows for reproducibility

Source: Story 6A.1 - Presence/Absence 矩阵生成, Story 6A.2 - 多物种比较
"""

import csv
import logging
from collections import defaultdict
from pathlib import Path

from workflow.lib.errors import CompGeneError, ErrorCode

logger = logging.getLogger(__name__)

# Module version for audit records
__version__ = "1.1.0"


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


# =============================================================================
# PA Matrix Reading (Story 6A.2)
# =============================================================================

def read_presence_absence_matrix(input_path: Path) -> tuple[dict[str, dict[str, int]], list[str]]:
    """
    Read presence/absence matrix TSV into dict and species list.

    Expected format (output of Story 6A.1 matrices_generate):
        orthogroup_id  species1  species2  species3
        OG0000000      1         1         0
        OG0000001      0         1         1

    Args:
        input_path: Path to presence_absence.tsv.

    Returns:
        Tuple of (matrix, species_list) where:
        - matrix: dict[orthogroup_id][species_id] = 0 or 1
        - species_list: ordered list of species from header

    Raises:
        CompGeneError: E_INPUT_MISSING if file does not exist.
        CompGeneError: E_INPUT_FORMAT if header is missing orthogroup_id column.
    """
    if not input_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Presence/absence matrix file not found: {input_path}",
        )

    matrix: dict[str, dict[str, int]] = {}

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None or "orthogroup_id" not in reader.fieldnames:
            raise CompGeneError(
                ErrorCode.E_INPUT_FORMAT,
                f"Invalid header in {input_path}. "
                f"Expected 'orthogroup_id' as first column, got: {reader.fieldnames}",
            )

        species_list = [col for col in reader.fieldnames if col != "orthogroup_id"]

        for row in reader:
            og_id = row["orthogroup_id"].strip()
            if not og_id:
                continue
            matrix[og_id] = {sp: int(row[sp]) for sp in species_list}

    logger.info(
        "Read PA matrix: %d orthogroups x %d species from %s",
        len(matrix),
        len(species_list),
        input_path.name,
    )
    return matrix, species_list


# =============================================================================
# Multi-Species Comparison (Story 6A.2)
# =============================================================================

def classify_orthogroup_sharing(
    matrix: dict[str, dict[str, int]],
    species_list: list[str],
) -> list[dict]:
    """
    Classify each orthogroup by sharing pattern across species.

    Categories:
    - all_shared: present in ALL species
    - species_specific: present in exactly 1 species
    - partial: present in >1 but <all species

    Args:
        matrix: dict[orthogroup_id][species_id] = 0 or 1.
        species_list: ordered list of species names.

    Returns:
        List of classification dicts sorted by orthogroup_id, each containing:
        - orthogroup_id: str
        - category: "all_shared" | "species_specific" | "partial"
        - n_species_present: int
        - present_species: comma-separated string
        - absent_species: comma-separated string
        - specific_to: species name if species_specific, else ""
    """
    n_total = len(species_list)
    classifications = []

    for og_id in sorted(matrix.keys()):
        row = matrix[og_id]
        present = [sp for sp in species_list if row.get(sp, 0) > 0]
        absent = [sp for sp in species_list if row.get(sp, 0) == 0]
        n_present = len(present)

        if n_present == 1:
            category = "species_specific"
        elif n_present == n_total:
            category = "all_shared"
        else:
            category = "partial"

        classifications.append({
            "orthogroup_id": og_id,
            "category": category,
            "n_species_present": n_present,
            "present_species": ",".join(present),
            "absent_species": ",".join(absent),
            "specific_to": present[0] if category == "species_specific" else "",
        })

    logger.info(
        "Classified %d orthogroups across %d species",
        len(classifications),
        n_total,
    )
    return classifications


def summarize_sharing_counts(
    classifications: list[dict],
    species_list: list[str],
) -> list[dict]:
    """
    Summarize orthogroup sharing classification counts.

    Args:
        classifications: Output of classify_orthogroup_sharing().
        species_list: Ordered list of species names.

    Returns:
        List of dicts with keys (metric, count), suitable for TSV output:
        - total_orthogroups
        - all_shared
        - partial_shared
        - species_specific_total
        - species_specific_{species_name} for each species
    """
    total = len(classifications)
    all_shared = sum(1 for c in classifications if c["category"] == "all_shared")
    partial = sum(1 for c in classifications if c["category"] == "partial")
    specific_total = sum(1 for c in classifications if c["category"] == "species_specific")

    per_species: dict[str, int] = {sp: 0 for sp in species_list}
    for c in classifications:
        if c["category"] == "species_specific":
            per_species[c["specific_to"]] += 1

    rows = [
        {"metric": "total_orthogroups", "count": total},
        {"metric": "all_shared", "count": all_shared},
        {"metric": "partial_shared", "count": partial},
        {"metric": "species_specific_total", "count": specific_total},
    ]
    for sp in species_list:
        rows.append({"metric": f"species_specific_{sp}", "count": per_species[sp]})

    return rows


def write_comparison_tsv(classifications: list[dict], output_path: Path) -> None:
    """
    Write orthogroup classification results to TSV.

    Uses atomic write (write to .tmp, then rename) for checkpoint safety.

    Args:
        classifications: Output of classify_orthogroup_sharing().
        output_path: Destination file path.
    """
    columns = [
        "orthogroup_id", "category", "n_species_present",
        "present_species", "absent_species", "specific_to",
    ]
    _atomic_write_tsv(classifications, columns, output_path)
    logger.info("Wrote %d classification rows to %s", len(classifications), output_path)


def write_summary_tsv(summary_rows: list[dict], output_path: Path) -> None:
    """
    Write comparison summary to TSV.

    Uses atomic write (write to .tmp, then rename) for checkpoint safety.

    Args:
        summary_rows: Output of summarize_sharing_counts().
        output_path: Destination file path.
    """
    columns = ["metric", "count"]
    _atomic_write_tsv(summary_rows, columns, output_path)
    logger.info("Wrote %d summary rows to %s", len(summary_rows), output_path)


def _atomic_write_tsv(rows: list[dict], columns: list[str], output_path: Path) -> None:
    """
    Write list of dicts to TSV with atomic write pattern.

    Args:
        rows: List of dicts to write.
        columns: Column names (keys in dicts) in output order.
        output_path: Destination file path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    try:
        with open(tmp_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for row in rows:
                writer.writerow([str(row.get(col, "")) for col in columns])
        tmp_path.rename(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
