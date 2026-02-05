"""
CompGene Annotation Enhancement Module.

This module provides functionality for enhancing lifted GFF3 annotations
with Liftoff classification attributes (coverage, identity, status).

Key features:
- Add liftoff_coverage, liftoff_identity, liftoff_status attributes
- Filter genes by classification status (present/uncertain/missing)
- Propagate attributes to child features (mRNA, CDS, exon)
- Atomic file writing for safety
- gzip compression support

Source: Story 5.3 - 迁移注释输出
"""

import csv
import logging
from pathlib import Path
from typing import Optional

from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.gff import open_gff_file, parse_gff3_attributes

logger = logging.getLogger(__name__)

# Module version for audit records
__version__ = "1.0.0"


# =============================================================================
# Constants
# =============================================================================

# Status values that are included in output by default
DEFAULT_ALLOWED_STATUSES = {"present"}

# Status values that are included when include_uncertain=True
UNCERTAIN_ALLOWED_STATUSES = {"present", "uncertain"}


# =============================================================================
# Classification Data Loading
# =============================================================================

def load_classification_data(classification_path: Path) -> dict[str, dict]:
    """
    Load gene classification data from TSV file.

    Args:
        classification_path: Path to gene_classification.tsv

    Returns:
        Dictionary mapping gene_id to classification data:
        {
            "gene1": {
                "status": "present",
                "coverage": 0.95,
                "identity": 0.92,
                "source_gene": "NM_001234"
            },
            ...
        }

    Raises:
        CompGeneError: If file not found or format invalid

    Example:
        >>> data = load_classification_data(Path("gene_classification.tsv"))
        >>> data["gene1"]["status"]
        'present'
    """
    if not classification_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Classification file not found: {classification_path}",
            details=str(classification_path)
        )

    result: dict[str, dict] = {}

    with open(classification_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            gene_id = row.get("gene_id", "").strip()
            if not gene_id:
                continue

            # Parse numeric values
            try:
                coverage = float(row.get("coverage", "0"))
            except (ValueError, TypeError):
                coverage = 0.0

            try:
                identity = float(row.get("identity", "0"))
            except (ValueError, TypeError):
                identity = 0.0

            # Normalize percentage values (>1.0) to 0.0-1.0 range
            if coverage > 1.0:
                coverage = coverage / 100.0
            if identity > 1.0:
                identity = identity / 100.0

            result[gene_id] = {
                "status": row.get("status", "").strip(),
                "coverage": coverage,
                "identity": identity,
                "source_gene": row.get("source_gene", "").strip(),
            }

    logger.info(f"Loaded classification data for {len(result)} genes")
    return result


# =============================================================================
# GFF3 Attribute Formatting
# =============================================================================

def format_gff3_attributes(attributes: dict[str, str]) -> str:
    """
    Format a dictionary of attributes into GFF3 format string.

    Args:
        attributes: Dictionary of attribute key-value pairs

    Returns:
        GFF3-formatted attribute string (key=value;key=value)

    Example:
        >>> attrs = {"ID": "gene1", "Name": "BRCA1"}
        >>> format_gff3_attributes(attrs)
        'ID=gene1;Name=BRCA1'
    """
    # URL-encode special characters
    def encode_value(v: str) -> str:
        return (v
                .replace(";", "%3B")
                .replace("=", "%3D")
                .replace("&", "%26")
                .replace(",", "%2C"))

    parts = []
    for key, value in attributes.items():
        if value:
            parts.append(f"{key}={encode_value(str(value))}")
        else:
            parts.append(f"{key}=")

    return ";".join(parts)


# =============================================================================
# GFF Line Enhancement
# =============================================================================

def enhance_gff_line(
    line: str,
    classification: dict,
    reference_species: str
) -> str:
    """
    Enhance a GFF3 line by adding liftoff classification attributes.

    Args:
        line: Original GFF3 line (tab-separated)
        classification: Classification data dict with status, coverage, identity, source_gene
        reference_species: Name of the reference species

    Returns:
        Enhanced GFF3 line with liftoff attributes added

    Example:
        >>> line = "chr1\\tliftoff\\tgene\\t1000\\t2000\\t.\\t+\\t.\\tID=gene1"
        >>> class_data = {"status": "present", "coverage": 0.95, "identity": 0.92, "source_gene": "NM_001234"}
        >>> enhanced = enhance_gff_line(line, class_data, "human")
        >>> "liftoff_status=present" in enhanced
        True
    """
    fields = line.strip().split("\t")
    if len(fields) != 9:
        return line

    # Parse existing attributes
    attrs = parse_gff3_attributes(fields[8])

    # Add liftoff classification attributes
    attrs["liftoff_coverage"] = str(classification.get("coverage", 0.0))
    attrs["liftoff_identity"] = str(classification.get("identity", 0.0))
    attrs["liftoff_status"] = classification.get("status", "")
    attrs["source_gene"] = classification.get("source_gene", "")
    attrs["reference_species"] = reference_species

    # Rebuild the line with enhanced attributes
    fields[8] = format_gff3_attributes(attrs)

    return "\t".join(fields)


# =============================================================================
# Status Filtering
# =============================================================================

def filter_by_status(
    statuses: list[str],
    include_uncertain: bool = False
) -> set[str]:
    """
    Determine which statuses should be included in output.

    Args:
        statuses: List of status values to filter
        include_uncertain: Whether to include 'uncertain' status

    Returns:
        Set of statuses that should be included

    Example:
        >>> filter_by_status(["present", "uncertain", "missing"], False)
        {'present'}
        >>> filter_by_status(["present", "uncertain", "missing"], True)
        {'present', 'uncertain'}
    """
    if include_uncertain:
        allowed = UNCERTAIN_ALLOWED_STATUSES
    else:
        allowed = DEFAULT_ALLOWED_STATUSES

    return {s for s in statuses if s in allowed}


# =============================================================================
# Main Enhancement Function
# =============================================================================

def enhance_gff_with_liftoff_attrs(
    gff_path: Path,
    classification_path: Path,
    output_path: Path,
    reference_species: str,
    include_uncertain: bool = False
) -> dict:
    """
    Enhance GFF3 file with Liftoff classification attributes.

    Reads the lifted annotation GFF3 and gene classification TSV,
    adds liftoff_* attributes to each feature, and filters to only
    include genes with the allowed status.

    Args:
        gff_path: Path to input lifted_annotation.gff3
        classification_path: Path to gene_classification.tsv
        output_path: Path for output enhanced GFF3
        reference_species: Name of the reference species
        include_uncertain: Include genes with 'uncertain' status (default: False)

    Returns:
        Statistics dictionary:
        - input_genes: Number of genes in input
        - output_genes: Number of genes in output
        - filtered_genes: Number of genes filtered out
        - unclassified_genes: Number of genes with no classification

    Raises:
        CompGeneError: If input files not found or invalid

    Example:
        >>> stats = enhance_gff_with_liftoff_attrs(
        ...     gff_path=Path("lifted.gff3"),
        ...     classification_path=Path("classification.tsv"),
        ...     output_path=Path("enhanced.gff3"),
        ...     reference_species="human"
        ... )
        >>> print(f"Output: {stats['output_genes']} genes")
    """
    # Validate input files
    if not gff_path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"GFF file not found: {gff_path}",
            details=str(gff_path)
        )

    # Load classification data
    classification_data = load_classification_data(classification_path)

    # Determine allowed statuses
    if include_uncertain:
        allowed_statuses = UNCERTAIN_ALLOWED_STATUSES
    else:
        allowed_statuses = DEFAULT_ALLOWED_STATUSES

    # Build gene_id -> classification mapping
    # Also build parent_id -> gene_id mapping for child features
    gene_classifications: dict[str, dict] = {}
    parent_to_gene: dict[str, str] = {}

    # Statistics
    input_genes = 0
    output_genes = 0
    filtered_genes = 0
    unclassified_genes = 0

    # First pass: identify which genes to include and build parent mapping
    genes_to_include: set[str] = set()

    with open_gff_file(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue

            feature_type = fields[2].lower()
            attrs = parse_gff3_attributes(fields[8])
            feature_id = attrs.get("ID", "")
            parent_id = attrs.get("Parent", "")

            if feature_type == "gene":
                input_genes += 1
                gene_id = feature_id

                if gene_id in classification_data:
                    class_info = classification_data[gene_id]
                    gene_classifications[gene_id] = class_info

                    if class_info["status"] in allowed_statuses:
                        genes_to_include.add(gene_id)
                        output_genes += 1
                    else:
                        filtered_genes += 1
                else:
                    # Gene not in classification - filter it out
                    unclassified_genes += 1

            elif parent_id:
                # Track parent-child relationships
                # mRNA's parent is gene, exon/CDS parent is mRNA
                parent_to_gene[feature_id] = parent_id

    # Build transitive parent mapping (child -> ultimate gene ancestor)
    def find_gene_ancestor(feature_id: str, visited: set) -> Optional[str]:
        if feature_id in visited:
            return None
        visited.add(feature_id)

        if feature_id in genes_to_include:
            return feature_id

        parent = parent_to_gene.get(feature_id)
        if parent:
            return find_gene_ancestor(parent, visited)
        return None

    # Second pass: write enhanced output
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open_gff_file(gff_path) as f_in, open(temp_path, "w") as f_out:
        for line in f_in:
            original_line = line
            line = line.strip()

            # Preserve header lines
            if not line or line.startswith("#"):
                f_out.write(original_line)
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue

            attrs = parse_gff3_attributes(fields[8])
            feature_id = attrs.get("ID", "")
            parent_id = attrs.get("Parent", "")
            feature_type = fields[2].lower()

            # Determine if this feature should be included
            gene_ancestor = None

            if feature_type == "gene":
                if feature_id in genes_to_include:
                    gene_ancestor = feature_id
            else:
                # Find the gene ancestor for child features
                gene_ancestor = find_gene_ancestor(parent_id, set())

            if gene_ancestor is None:
                # Skip features not associated with included genes
                continue

            # Get classification for the gene ancestor
            class_info = gene_classifications.get(gene_ancestor)
            if class_info is None:
                continue

            # Enhance the line
            enhanced_line = enhance_gff_line(line, class_info, reference_species)
            f_out.write(enhanced_line + "\n")

    # Atomic rename
    temp_path.rename(output_path)

    stats = {
        "input_genes": input_genes,
        "output_genes": output_genes,
        "filtered_genes": filtered_genes,
        "unclassified_genes": unclassified_genes,
    }

    logger.info(f"Enhanced GFF: {input_genes} input genes, {output_genes} output, "
                f"{filtered_genes} filtered, {unclassified_genes} unclassified")

    return stats
