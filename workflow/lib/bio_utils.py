"""
CompGene Biopython Utilities Module.

This module provides Biopython-enhanced sequence operations with graceful
degradation when Biopython is not available.

Key features:
- Multiple genetic code table support for translation
- Multi-format FASTA ID parsing (UniProt, NCBI, OrthoFinder, etc.)
- ID consistency checking between pipeline stages
- Audit summary generation (summary.json)

Source: Story 2.2b - Biopython 工具封装
"""

import json
import re
import statistics
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.fasta import parse_fasta, get_sequence_stats, SequenceStats


# =============================================================================
# Biopython Availability Check
# =============================================================================

_BIOPYTHON_AVAILABLE: Optional[bool] = None


def check_biopython_available() -> bool:
    """
    Check if Biopython is available for import.

    Returns:
        True if Biopython can be imported, False otherwise.

    Note:
        Result is cached after first check for performance.
    """
    global _BIOPYTHON_AVAILABLE
    if _BIOPYTHON_AVAILABLE is None:
        try:
            from Bio import SeqIO  # noqa: F401
            from Bio.Seq import Seq  # noqa: F401
            _BIOPYTHON_AVAILABLE = True
        except ImportError:
            _BIOPYTHON_AVAILABLE = False
    return _BIOPYTHON_AVAILABLE


def require_biopython(func_name: str) -> None:
    """
    Raise error if Biopython is not available.

    Args:
        func_name: Name of the function requiring Biopython.

    Raises:
        CompGeneError: If Biopython is not installed.
    """
    if not check_biopython_available():
        raise CompGeneError(
            ErrorCode.E_TOOL_NOT_FOUND,
            f"Biopython is required for {func_name}",
            details="Install with: pip install biopython>=1.81"
        )


# =============================================================================
# Genetic Code Tables
# =============================================================================

CODON_TABLE_NAMES: dict[int, str] = {
    1: "Standard",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold Mitochondrial",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate Nuclear",
    9: "Echinoderm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Alternative Flatworm Mitochondrial",
    15: "Blepharisma Macronuclear",
    16: "Chlorophycean Mitochondrial",
    21: "Trematode Mitochondrial",
    22: "Scenedesmus obliquus Mitochondrial",
    23: "Thraustochytrium Mitochondrial",
}

# Standard genetic code for fallback (table_id=1)
STANDARD_CODON_TABLE: dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_with_table(
    sequence: str,
    table_id: int = 1,
    to_stop: bool = False,
    stop_symbol: str = "*"
) -> str:
    """
    Translate a nucleotide sequence using specified genetic code table.

    Args:
        sequence: DNA/RNA sequence to translate.
        table_id: NCBI genetic code table ID (default: 1 = Standard).
        to_stop: If True, stop translation at first stop codon.
        stop_symbol: Symbol to use for stop codons (default: "*").

    Returns:
        Translated amino acid sequence.

    Raises:
        CompGeneError: If table_id is invalid or Biopython required but unavailable.

    Note:
        For table_id=1 (Standard), falls back to built-in table if Biopython
        is not available. Other tables require Biopython.
    """
    if table_id not in CODON_TABLE_NAMES:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"Invalid genetic code table ID: {table_id}",
            details=f"Valid table IDs: {list(CODON_TABLE_NAMES.keys())}"
        )

    sequence = sequence.upper().replace('U', 'T')

    # Use Biopython if available
    if check_biopython_available():
        from Bio.Seq import Seq
        seq = Seq(sequence)
        protein = str(seq.translate(table=table_id, to_stop=to_stop, stop_symbol=stop_symbol))
        return protein

    # Fallback for standard table only
    if table_id != 1:
        require_biopython(f"translate_with_table(table_id={table_id})")

    # Manual translation with standard table
    protein_parts = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if 'N' in codon or len(codon) < 3:
            aa = 'X'
        else:
            aa = STANDARD_CODON_TABLE.get(codon, 'X')

        if aa == '*':
            if to_stop:
                break
            protein_parts.append(stop_symbol)
        else:
            protein_parts.append(aa)

    return ''.join(protein_parts)


# =============================================================================
# ID Parsing
# =============================================================================

from enum import Enum


class IDFormat(str, Enum):
    """Enumeration for FASTA ID format detection."""
    UNIPROT = "uniprot"          # sp|P12345|GENE_HUMAN or tr|A0A123|...
    NCBI_GI = "ncbi_gi"          # gi|123456|ref|NP_001234.1|
    NCBI_REF = "ncbi_ref"        # ref|NP_001234.1|
    GENBANK = "genbank"          # gb|AAA12345.1|
    ORTHOFINDER = "orthofinder"  # species|gene_id
    SIMPLE = "simple"            # gene_id or gene_id description
    UNKNOWN = "unknown"


def detect_id_format(header: str) -> IDFormat:
    """
    Detect the format of a FASTA header.

    Args:
        header: FASTA header line (with or without ">").

    Returns:
        One of IDFormat constants indicating detected format.
    """
    header = header.lstrip(">").strip()

    if not header:
        return IDFormat.UNKNOWN

    # UniProt: sp|P12345|GENE_HUMAN or tr|A0A123|...
    if header.startswith(("sp|", "tr|")):
        return IDFormat.UNIPROT

    # NCBI GI: gi|123456|ref|NP_001234.1|
    if header.startswith("gi|"):
        return IDFormat.NCBI_GI

    # NCBI RefSeq: ref|NP_001234.1|
    if header.startswith("ref|"):
        return IDFormat.NCBI_REF

    # GenBank: gb|AAA12345.1|
    if header.startswith("gb|"):
        return IDFormat.GENBANK

    # OrthoFinder style: species|gene_id (exactly one pipe, no spaces before pipe)
    if "|" in header:
        parts = header.split("|")
        if len(parts) == 2 and " " not in parts[0]:
            return IDFormat.ORTHOFINDER

    # Simple format: gene_id or gene_id description
    return IDFormat.SIMPLE


def parse_fasta_id(header: str, format_hint: Optional[str] = None) -> str:
    """
    Parse a FASTA header to extract the primary sequence ID.

    Args:
        header: FASTA header line (with or without ">").
        format_hint: Optional format hint (one of IDFormat constants).
            If None, format is auto-detected.

    Returns:
        Extracted sequence ID.

    Examples:
        >>> parse_fasta_id("sp|P12345|GENE_HUMAN")
        'P12345'
        >>> parse_fasta_id("gi|123|ref|NP_001234.1|")
        'NP_001234.1'
        >>> parse_fasta_id("mmur|MMUR_00001")
        'MMUR_00001'
        >>> parse_fasta_id("gene_001 some description")
        'gene_001'
    """
    header = header.lstrip(">").strip()

    if not header:
        return ""

    fmt = format_hint or detect_id_format(header)

    if fmt == IDFormat.UNIPROT:
        # sp|P12345|GENE_HUMAN -> P12345
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1]

    elif fmt == IDFormat.NCBI_GI:
        # gi|123456|ref|NP_001234.1| -> NP_001234.1
        parts = header.split("|")
        # Find "ref" and take the next part
        for i, part in enumerate(parts):
            if part == "ref" and i + 1 < len(parts):
                return parts[i + 1].rstrip("|")
        # Fallback: take second part
        if len(parts) >= 2:
            return parts[1]

    elif fmt == IDFormat.NCBI_REF:
        # ref|NP_001234.1| -> NP_001234.1
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1].rstrip("|")

    elif fmt == IDFormat.GENBANK:
        # gb|AAA12345.1| -> AAA12345.1
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1].rstrip("|")

    elif fmt == IDFormat.ORTHOFINDER:
        # species|gene_id -> gene_id
        parts = header.split("|")
        if len(parts) >= 2:
            # Take last part, strip description
            last_part = parts[-1]
            return last_part.split()[0] if " " in last_part else last_part

    # Simple or unknown: take first whitespace-delimited token
    return header.split()[0]


def extract_ids(
    fasta_path: Path,
    format_hint: Optional[str] = None,
    include_raw: bool = False
) -> set[str]:
    """
    Extract all sequence IDs from a FASTA file.

    Args:
        fasta_path: Path to FASTA file (supports .gz).
        format_hint: Optional format hint for ID parsing.
        include_raw: If True, also include raw IDs (first token).

    Returns:
        Set of extracted sequence IDs.

    Raises:
        CompGeneError: If file not found or invalid format.
    """
    ids: set[str] = set()

    for record in parse_fasta(fasta_path):
        parsed_id = parse_fasta_id(record.description, format_hint)
        ids.add(parsed_id)

        if include_raw:
            ids.add(record.id)

    return ids


def extract_ids_with_mapping(
    fasta_path: Path,
    format_hint: Optional[str] = None
) -> tuple[set[str], dict[str, str]]:
    """
    Extract IDs and create a mapping from raw ID to parsed ID.

    Args:
        fasta_path: Path to FASTA file.
        format_hint: Optional format hint for ID parsing.

    Returns:
        Tuple of (set of parsed IDs, dict mapping raw_id -> parsed_id).
    """
    ids: set[str] = set()
    mapping: dict[str, str] = {}

    for record in parse_fasta(fasta_path):
        parsed_id = parse_fasta_id(record.description, format_hint)
        ids.add(parsed_id)
        mapping[record.id] = parsed_id

    return ids, mapping


# =============================================================================
# ID Consistency Checking
# =============================================================================

@dataclass
class ComparisonResult:
    """Result of comparing two ID sets."""

    source_a: str
    source_b: str
    total_a: int
    total_b: int
    common: int
    only_in_a: set[str] = field(default_factory=set)
    only_in_b: set[str] = field(default_factory=set)

    @property
    def match_rate(self) -> float:
        """Calculate match rate (common / total_a)."""
        if self.total_a == 0:
            return 0.0
        return self.common / self.total_a

    def is_consistent(self, threshold: float = 0.95) -> bool:
        """Check if match rate meets threshold."""
        return self.match_rate >= threshold

    def to_dict(self, max_items: int = 100) -> dict:
        """Convert to dictionary for JSON serialization."""
        only_in_a_list = sorted(list(self.only_in_a))
        only_in_b_list = sorted(list(self.only_in_b))

        return {
            "source_a": self.source_a,
            "source_b": self.source_b,
            "comparison": {
                "total_a": self.total_a,
                "total_b": self.total_b,
                "common": self.common,
                "only_in_a_count": len(self.only_in_a),
                "only_in_b_count": len(self.only_in_b),
                "match_rate": round(self.match_rate, 6),
            },
            "only_in_a": only_in_a_list[:max_items],
            "only_in_a_truncated": len(only_in_a_list) > max_items,
            "only_in_b": only_in_b_list[:max_items],
            "only_in_b_truncated": len(only_in_b_list) > max_items,
        }


def compare_id_sets(
    ids_a: set[str],
    ids_b: set[str],
    source_a: str = "source_a",
    source_b: str = "source_b"
) -> ComparisonResult:
    """
    Compare two sets of IDs.

    Args:
        ids_a: First set of IDs.
        ids_b: Second set of IDs.
        source_a: Description of first source.
        source_b: Description of second source.

    Returns:
        ComparisonResult with comparison statistics.
    """
    common = ids_a & ids_b
    only_in_a = ids_a - ids_b
    only_in_b = ids_b - ids_a

    return ComparisonResult(
        source_a=source_a,
        source_b=source_b,
        total_a=len(ids_a),
        total_b=len(ids_b),
        common=len(common),
        only_in_a=only_in_a,
        only_in_b=only_in_b,
    )


def extract_ids_from_tsv(
    tsv_path: Path,
    id_column: int = 0,
    skip_header: bool = True,
    delimiter: str = "\t"
) -> set[str]:
    """
    Extract IDs from a TSV file column.

    Args:
        tsv_path: Path to TSV file.
        id_column: Column index containing IDs (0-based).
        skip_header: Whether to skip the first line.
        delimiter: Column delimiter.

    Returns:
        Set of IDs from the specified column.
    """
    import gzip

    ids: set[str] = set()

    open_func = gzip.open if str(tsv_path).endswith('.gz') else open
    with open_func(tsv_path, 'rt', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if skip_header and i == 0:
                continue
            line = line.strip()
            if not line:
                continue
            parts = line.split(delimiter)
            if id_column < len(parts):
                cell = parts[id_column].strip()
                # Handle comma-separated IDs in cell
                if ',' in cell:
                    for sub_id in cell.split(','):
                        sub_id = sub_id.strip()
                        if sub_id:
                            ids.add(sub_id)
                elif cell:
                    ids.add(cell)

    return ids


def check_orthofinder_consistency(
    input_fasta: Path,
    orthogroups_tsv: Path,
    species_id: Optional[str] = None
) -> ComparisonResult:
    """
    Check ID consistency between input FASTA and OrthoFinder output.

    Args:
        input_fasta: Path to input protein FASTA.
        orthogroups_tsv: Path to OrthoFinder Orthogroups.tsv.
        species_id: Species ID to filter (if None, extracts all).

    Returns:
        ComparisonResult comparing input IDs to OrthoFinder output IDs.
    """
    # Extract IDs from input FASTA
    input_ids = extract_ids(input_fasta)

    # Extract IDs from OrthoFinder output
    # Orthogroups.tsv format: OG_id\tspecies1_genes\tspecies2_genes\t...
    orthofinder_ids: set[str] = set()

    import gzip
    open_func = gzip.open if str(orthogroups_tsv).endswith('.gz') else open

    with open_func(orthogroups_tsv, 'rt', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        species_columns = header[1:]  # Skip Orthogroup column

        # Find column index for species if specified
        col_indices = list(range(1, len(header)))
        if species_id:
            col_indices = [i for i, s in enumerate(header) if s == species_id]

        for line in f:
            parts = line.strip().split('\t')
            for col_idx in col_indices:
                if col_idx < len(parts):
                    cell = parts[col_idx].strip()
                    if cell:
                        # OrthoFinder uses comma-separated gene IDs
                        for gene_id in cell.split(','):
                            gene_id = gene_id.strip()
                            if gene_id:
                                # Parse the ID to get core part
                                orthofinder_ids.add(parse_fasta_id(gene_id))

    return compare_id_sets(
        input_ids,
        orthofinder_ids,
        source_a=str(input_fasta),
        source_b=str(orthogroups_tsv)
    )


def check_eggnog_consistency(
    input_fasta: Path,
    eggnog_annotations: Path
) -> ComparisonResult:
    """
    Check ID consistency between input FASTA and eggNOG-mapper output.

    Args:
        input_fasta: Path to input protein FASTA.
        eggnog_annotations: Path to eggNOG-mapper annotations file.

    Returns:
        ComparisonResult comparing input IDs to eggNOG annotated IDs.
    """
    # Extract IDs from input FASTA
    input_ids = extract_ids(input_fasta)

    # Extract IDs from eggNOG output (first column, skip comment lines)
    eggnog_ids: set[str] = set()

    import gzip
    open_func = gzip.open if str(eggnog_annotations).endswith('.gz') else open

    with open_func(eggnog_annotations, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts:
                gene_id = parts[0].strip()
                if gene_id:
                    eggnog_ids.add(parse_fasta_id(gene_id))

    return compare_id_sets(
        input_ids,
        eggnog_ids,
        source_a=str(input_fasta),
        source_b=str(eggnog_annotations)
    )


def write_consistency_report(
    result: ComparisonResult,
    output_path: Path,
    threshold: float = 0.95
) -> None:
    """
    Write consistency check result to JSON file.

    Args:
        result: ComparisonResult to write.
        output_path: Path for output JSON file.
        threshold: Threshold for PASS/FAIL status.
    """
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        **result.to_dict(),
        "status": "PASS" if result.is_consistent(threshold) else "FAIL",
        "threshold": threshold,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


# =============================================================================
# Audit Summary Generation
# =============================================================================

def compute_length_histogram(
    lengths: list[int],
    bins: Optional[list[int]] = None
) -> dict[str, list]:
    """
    Compute length distribution histogram.

    Args:
        lengths: List of sequence lengths.
        bins: Bin edges (default: [0, 100, 200, 500, 1000, 2000, 5000, 10000, inf]).

    Returns:
        Dict with 'bins' and 'counts' lists.
    """
    if bins is None:
        bins = [0, 100, 200, 500, 1000, 2000, 5000, 10000]

    # Add infinity as last bin edge
    bin_edges = bins + [float('inf')]
    counts = [0] * (len(bin_edges) - 1)

    for length in lengths:
        for i in range(len(bin_edges) - 1):
            if bin_edges[i] <= length < bin_edges[i + 1]:
                counts[i] += 1
                break

    return {
        "bins": bins,
        "counts": counts
    }


def detect_sequence_type(sequences: list[str], sample_size: int = 100) -> str:
    """
    Detect whether sequences are nucleotide or protein.

    Args:
        sequences: List of sequences to analyze.
        sample_size: Number of sequences to sample.

    Returns:
        "nucleotide" or "protein".
    """
    if not sequences:
        return "unknown"

    # Sample sequences
    sample = sequences[:sample_size]
    combined = ''.join(sample).upper()

    if not combined:
        return "unknown"

    # Count nucleotide bases
    nuc_bases = set("ATGCUN")
    nuc_count = sum(1 for c in combined if c in nuc_bases)
    nuc_ratio = nuc_count / len(combined)

    return "nucleotide" if nuc_ratio > 0.9 else "protein"


def generate_audit_summary(fasta_path: Path) -> dict:
    """
    Generate audit summary statistics for a FASTA file.

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        Dictionary with comprehensive statistics.
    """
    # Collect data
    lengths: list[int] = []
    sequences: list[str] = []
    ids: list[str] = []
    gc_total = 0
    n_total = 0
    x_total = 0  # For proteins

    for record in parse_fasta(fasta_path):
        lengths.append(len(record.sequence))
        sequences.append(record.sequence)
        ids.append(record.id)

        seq_upper = record.sequence.upper()
        gc_total += seq_upper.count('G') + seq_upper.count('C')
        n_total += seq_upper.count('N')
        x_total += seq_upper.count('X')

    if not lengths:
        return {
            "file_path": fasta_path.name,
            "file_size_bytes": fasta_path.stat().st_size if fasta_path.exists() else 0,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "sequence_type": "unknown",
            "status": "empty",
            "warning": "FASTA file contains no sequences",
            "statistics": {
                "sequence_count": 0,
                "total_length": 0,
            }
        }

    total_length = sum(lengths)
    seq_type = detect_sequence_type(sequences)

    # Calculate GC content for nucleotides
    gc_content = None
    n_content = 0.0
    if seq_type == "nucleotide" and total_length > 0:
        atgc_total = sum(
            s.upper().count('A') + s.upper().count('T') +
            s.upper().count('G') + s.upper().count('C')
            for s in sequences
        )
        if atgc_total > 0:
            gc_content = round(gc_total / atgc_total, 6)
        n_content = round(n_total / total_length, 6)
    elif seq_type == "protein" and total_length > 0:
        n_content = round(x_total / total_length, 6)  # X is ambiguous for proteins

    # Calculate N50
    sorted_lengths = sorted(lengths, reverse=True)
    half_total = total_length / 2
    cumsum = 0
    n50 = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= half_total:
            n50 = length
            break

    # Detect ID format
    sample_ids = ids[:3]
    detected_format = detect_id_format(ids[0]) if ids else IDFormat.UNKNOWN

    return {
        "file_path": fasta_path.name,
        "file_size_bytes": fasta_path.stat().st_size if fasta_path.exists() else 0,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sequence_type": seq_type,
        "statistics": {
            "sequence_count": len(lengths),
            "total_length": total_length,
            "mean_length": round(total_length / len(lengths), 2),
            "median_length": int(statistics.median(lengths)),
            "min_length": min(lengths),
            "max_length": max(lengths),
            "n50": n50,
            "gc_content": gc_content,
            "n_content": n_content,
            "ambiguous_count": n_total if seq_type == "nucleotide" else x_total,
            "length_histogram": compute_length_histogram(lengths),
        },
        "id_format": {
            "detected": detected_format,
            "sample_ids": sample_ids,
        }
    }


def write_summary_json(
    fasta_path: Path,
    output_path: Optional[Path] = None
) -> Path:
    """
    Generate and write audit summary to JSON file.

    Args:
        fasta_path: Path to FASTA file.
        output_path: Path for output JSON. If None, uses fasta_path with
            .summary.json extension.

    Returns:
        Path to written summary file.
    """
    if output_path is None:
        # Replace .fa.gz or .fasta.gz with .summary.json
        name = fasta_path.name
        if name.endswith('.fa.gz'):
            name = name[:-6] + '.summary.json'
        elif name.endswith('.fasta.gz'):
            name = name[:-9] + '.summary.json'
        elif name.endswith('.fa'):
            name = name[:-3] + '.summary.json'
        elif name.endswith('.fasta'):
            name = name[:-6] + '.summary.json'
        else:
            name = name + '.summary.json'
        output_path = fasta_path.parent / name

    summary = generate_audit_summary(fasta_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    return output_path
