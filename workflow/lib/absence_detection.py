"""
CompGene Absence Detection Module.

This module provides functionality for analyzing Liftoff mapping quality
and classifying genes as present, uncertain, or missing based on
coverage and identity thresholds.

Key features:
- Gene classification based on configurable thresholds
- Coverage/identity extraction from Liftoff GFF attributes
- Support for unmapped genes from Liftoff output
- Statistics calculation for classification results
- gzip compression support

Source: Story 5.2 - 映射质量分析与缺失检测
"""

import logging
from pathlib import Path
from typing import Optional

from workflow.lib.errors import ErrorCode, CompGeneError

logger = logging.getLogger(__name__)
from workflow.lib.gff import parse_gff3

# Module version for audit records
__version__ = "1.0.0"


# =============================================================================
# Constants
# =============================================================================

DEFAULT_MIN_COVERAGE = 0.50
DEFAULT_MIN_IDENTITY = 0.50


# =============================================================================
# Gene Classification
# =============================================================================

def classify_gene(
    coverage: float,
    identity: Optional[float],
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    min_identity: float = DEFAULT_MIN_IDENTITY
) -> str:
    """
    Classify a gene's status based on coverage and identity.

    Classification rules:
    - present: coverage >= min_coverage AND identity >= min_identity
    - uncertain: coverage > 0 but doesn't meet thresholds
    - missing: coverage == 0

    Args:
        coverage: Liftoff coverage value (0.0-1.0)
        identity: Sequence identity value (0.0-1.0), defaults to coverage if None
        min_coverage: Minimum coverage threshold (default: 0.50)
        min_identity: Minimum identity threshold (default: 0.50)

    Returns:
        Classification string: 'present', 'uncertain', or 'missing'

    Example:
        >>> classify_gene(coverage=0.95, identity=0.92)
        'present'
        >>> classify_gene(coverage=0.42, identity=0.38)
        'uncertain'
        >>> classify_gene(coverage=0.0, identity=0.0)
        'missing'
    """
    # Default identity to coverage if not provided
    if identity is None:
        identity = coverage

    # Classification logic
    if coverage >= min_coverage and identity >= min_identity:
        return "present"
    elif coverage > 0:
        return "uncertain"
    else:
        return "missing"


# =============================================================================
# Coverage/Identity Extraction
# =============================================================================

def extract_coverage_identity(attributes: dict[str, str]) -> tuple[float, float]:
    """
    Extract coverage and identity values from GFF attribute dictionary.

    Handles various attribute naming conventions:
    - coverage: Liftoff standard attribute
    - identity / sequence_identity: Identity attributes

    If identity is not present, it defaults to the coverage value.
    If values > 1.0, they are assumed to be percentages and normalized.

    Args:
        attributes: Dictionary of GFF attribute key-value pairs

    Returns:
        Tuple of (coverage, identity) as floats (0.0-1.0)

    Example:
        >>> attrs = {"coverage": "0.95", "identity": "0.92"}
        >>> coverage, identity = extract_coverage_identity(attrs)
        >>> print(f"coverage={coverage}, identity={identity}")
        coverage=0.95, identity=0.92
    """
    # Extract coverage
    coverage_str = attributes.get("coverage", "0")
    try:
        coverage = float(coverage_str)
    except (ValueError, TypeError):
        coverage = 0.0

    # Normalize if percentage (> 1.0)
    if coverage > 1.0:
        coverage = coverage / 100.0

    # Extract identity (try multiple attribute names)
    # Note: Liftoff uses "sequence_ID" for sequence identity score
    identity = None
    for key in ("identity", "sequence_identity", "sequence_ID"):
        if key in attributes:
            try:
                identity = float(attributes[key])
                # Normalize if percentage
                if identity > 1.0:
                    identity = identity / 100.0
                break
            except (ValueError, TypeError):
                pass

    # Default identity to coverage if not found
    if identity is None:
        identity = coverage

    return coverage, identity


# =============================================================================
# Liftoff Gene Classification
# =============================================================================

def classify_liftoff_genes(
    lifted_gff: Path,
    unmapped_file: Path,
    reference_species: str,
    target_species: str,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    min_identity: float = DEFAULT_MIN_IDENTITY
) -> list[dict]:
    """
    Classify genes from Liftoff output files.

    Parses the lifted GFF annotation and unmapped features file to
    classify each gene as present, uncertain, or missing.

    Args:
        lifted_gff: Path to Liftoff lifted_annotation.gff3
        unmapped_file: Path to Liftoff unmapped_features.txt
        reference_species: Name of the reference species
        target_species: Name of the target species
        min_coverage: Minimum coverage threshold (default: 0.50)
        min_identity: Minimum identity threshold (default: 0.50)

    Returns:
        List of dictionaries with classification results.
        Each dict contains:
        - gene_id: Gene identifier
        - gene_name: Gene name (if available)
        - reference_species: Reference species name
        - target_species: Target species name
        - status: Classification ('present', 'uncertain', 'missing')
        - coverage: Coverage value
        - identity: Identity value
        - source_gene: Original reference gene ID

    Raises:
        CompGeneError: E_INPUT_MISSING if files don't exist

    Example:
        >>> results = classify_liftoff_genes(
        ...     lifted_gff=Path("lifted.gff3"),
        ...     unmapped_file=Path("unmapped.txt"),
        ...     reference_species="human",
        ...     target_species="mouse_lemur"
        ... )
    """
    # Validate inputs
    if not lifted_gff.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Lifted GFF file not found: {lifted_gff}",
            details=str(lifted_gff)
        )

    results = []

    # Parse lifted GFF - only gene features
    for feature in parse_gff3(lifted_gff):
        if feature.type.lower() != "gene":
            continue

        coverage, identity = extract_coverage_identity(feature.attributes)
        status = classify_gene(coverage, identity, min_coverage, min_identity)

        results.append({
            "gene_id": feature.id or feature.attributes.get("ID", ""),
            "gene_name": feature.name or feature.attributes.get("Name", ""),
            "reference_species": reference_species,
            "target_species": target_species,
            "status": status,
            "coverage": coverage,
            "identity": identity,
            "source_gene": feature.id or feature.attributes.get("ID", ""),
        })

    # Parse unmapped features - all are classified as missing
    if not unmapped_file.exists():
        logger.warning(f"Unmapped features file not found: {unmapped_file}. "
                      "Assuming no unmapped genes.")
    elif unmapped_file.exists():
        with open(unmapped_file, "r") as f:
            for line in f:
                gene_id = line.strip()
                if not gene_id:
                    continue

                results.append({
                    "gene_id": gene_id,
                    "gene_name": "",
                    "reference_species": reference_species,
                    "target_species": target_species,
                    "status": "missing",
                    "coverage": 0.0,
                    "identity": 0.0,
                    "source_gene": gene_id,
                })

    return results


# =============================================================================
# Statistics Calculation
# =============================================================================

def calculate_classification_stats(results: list[dict]) -> dict:
    """
    Calculate summary statistics from classification results.

    Args:
        results: List of classification result dictionaries

    Returns:
        Dictionary with statistics:
        - total_genes: Total number of genes
        - present: Count of present genes
        - uncertain: Count of uncertain genes
        - missing: Count of missing genes
        - present_rate: Fraction of present genes
        - missing_rate: Fraction of missing genes

    Example:
        >>> stats = calculate_classification_stats(results)
        >>> print(f"Present rate: {stats['present_rate']:.1%}")
    """
    total = len(results)

    if total == 0:
        return {
            "total_genes": 0,
            "present": 0,
            "uncertain": 0,
            "missing": 0,
            "present_rate": 0.0,
            "missing_rate": 0.0,
        }

    present = sum(1 for r in results if r["status"] == "present")
    uncertain = sum(1 for r in results if r["status"] == "uncertain")
    missing = sum(1 for r in results if r["status"] == "missing")

    return {
        "total_genes": total,
        "present": present,
        "uncertain": uncertain,
        "missing": missing,
        "present_rate": round(present / total, 4),
        "missing_rate": round(missing / total, 4),
    }


# =============================================================================
# Filtering Functions
# =============================================================================

def filter_missing_genes(
    results: list[dict],
    include_uncertain: bool = False
) -> list[dict]:
    """
    Filter classification results to only missing genes.

    Args:
        results: List of classification result dictionaries
        include_uncertain: If True, include uncertain genes in output

    Returns:
        Filtered list of missing (and optionally uncertain) genes

    Example:
        >>> missing = filter_missing_genes(results)
        >>> print(f"Found {len(missing)} missing genes")
    """
    if include_uncertain:
        return [r for r in results if r["status"] in ("missing", "uncertain")]
    else:
        return [r for r in results if r["status"] == "missing"]


def read_classification_stats_from_tsv(classification_path: Path) -> dict:
    """
    Read classification TSV and calculate statistics.

    This is a convenience function for summarizing already-written
    classification results without re-running classification.

    Args:
        classification_path: Path to gene_classification.tsv

    Returns:
        Dictionary with statistics (same format as calculate_classification_stats)

    Example:
        >>> stats = read_classification_stats_from_tsv(Path("gene_classification.tsv"))
        >>> print(f"Present: {stats['present']}, Missing: {stats['missing']}")
    """
    import csv

    counts = {"present": 0, "uncertain": 0, "missing": 0}

    if not classification_path.exists():
        return {
            "total_genes": 0,
            "present": 0,
            "uncertain": 0,
            "missing": 0,
            "present_rate": 0.0,
            "missing_rate": 0.0,
        }

    with open(classification_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            status = row.get("status", "")
            if status in counts:
                counts[status] += 1

    total = sum(counts.values())

    return {
        "total_genes": total,
        "present": counts["present"],
        "uncertain": counts["uncertain"],
        "missing": counts["missing"],
        "present_rate": round(counts["present"] / total, 4) if total > 0 else 0.0,
        "missing_rate": round(counts["missing"] / total, 4) if total > 0 else 0.0,
    }
