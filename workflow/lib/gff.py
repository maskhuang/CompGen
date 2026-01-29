"""
CompGene GFF3/GTF Parser Module.

This module provides streaming parsers for GFF3 and GTF annotation formats,
with support for building gene hierarchies and selecting representative transcripts.

Key features:
- Streaming/iterator-based parsing for memory efficiency (NFR2)
- GFF3 and GTF format support with auto-detection
- Gene → Transcript → CDS hierarchy construction
- Representative transcript selection (longest CDS/exon)
- gzip compression support
- Format validation with detailed error reporting

Source: Story 2.1 - GFF/FASTA 解析库
"""

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Iterator, Optional

from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class GFFFeature:
    """
    Represents a single feature from a GFF3/GTF file.

    Corresponds to one line in the annotation file.

    Attributes:
        seqid: Sequence identifier (chromosome/contig name).
        source: Feature source (e.g., "NCBI", "ensembl").
        type: Feature type (e.g., "gene", "mRNA", "CDS", "exon").
        start: Start position (1-based, inclusive).
        end: End position (1-based, inclusive).
        score: Feature score (float or None if ".").
        strand: Strand ("+" or "-" or "." for unstranded).
        phase: CDS phase (0, 1, 2, or None for non-CDS features).
        attributes: Dictionary of attribute key-value pairs.
        line_num: Line number in source file (for error reporting).
    """

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: Optional[float]
    strand: str
    phase: Optional[int]
    attributes: dict[str, str]
    line_num: int = 0

    @property
    def id(self) -> Optional[str]:
        """Get the feature ID (GFF3) or construct from gene_id/transcript_id (GTF)."""
        return self.attributes.get("ID") or self.attributes.get("transcript_id")

    @property
    def parent(self) -> Optional[str]:
        """Get the Parent attribute (GFF3) or gene_id (GTF)."""
        return self.attributes.get("Parent") or self.attributes.get("gene_id")

    @property
    def name(self) -> Optional[str]:
        """Get the feature Name attribute."""
        return self.attributes.get("Name") or self.attributes.get("gene_name")

    @property
    def length(self) -> int:
        """Get the feature length in base pairs."""
        return self.end - self.start + 1


@dataclass
class TranscriptModel:
    """
    Represents a transcript with its exons and CDS features.

    Attributes:
        transcript_id: Unique transcript identifier.
        gene_id: Parent gene identifier.
        exons: List of exon features.
        cds_features: List of CDS features.
        feature: The original mRNA/transcript feature.
    """

    transcript_id: str
    gene_id: str
    exons: list[GFFFeature] = field(default_factory=list)
    cds_features: list[GFFFeature] = field(default_factory=list)
    feature: Optional[GFFFeature] = None

    @property
    def cds_length(self) -> int:
        """Total CDS length in base pairs."""
        return sum(cds.length for cds in self.cds_features)

    @property
    def exon_length(self) -> int:
        """Total exon length in base pairs."""
        return sum(exon.length for exon in self.exons)

    @property
    def effective_length(self) -> int:
        """Effective length for comparison (CDS if available, else exon)."""
        return self.cds_length if self.cds_features else self.exon_length


@dataclass
class GeneModel:
    """
    Represents a gene with its transcripts and hierarchy.

    Attributes:
        gene_id: Unique gene identifier.
        transcripts: Dictionary of transcript_id -> TranscriptModel.
        representative_transcript: ID of the representative transcript (longest).
        feature: The original gene feature.
    """

    gene_id: str
    transcripts: dict[str, TranscriptModel] = field(default_factory=dict)
    representative_transcript: Optional[str] = None
    feature: Optional[GFFFeature] = None

    def get_representative(self) -> Optional[TranscriptModel]:
        """Get the representative transcript model."""
        if self.representative_transcript:
            return self.transcripts.get(self.representative_transcript)
        return None


# =============================================================================
# File Utilities
# =============================================================================

def open_gff_file(path: Path, mode: str = "rt") -> IO[str]:
    """
    Open a GFF/GTF file, automatically detecting gzip compression.

    Args:
        path: Path to the file.
        mode: File open mode (default: "rt" for text reading).

    Returns:
        File handle for reading.

    Raises:
        CompGeneError: If file cannot be opened.
    """
    if not path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"GFF/GTF file not found: {path}",
            details=str(path)
        )

    try:
        if path.suffix == ".gz" or str(path).endswith(".gz"):
            return gzip.open(path, mode, encoding="utf-8")
        return open(path, mode, encoding="utf-8")
    except Exception as e:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"Cannot open file: {path}",
            details=str(e)
        )


# =============================================================================
# GFF3 Parsing
# =============================================================================

def parse_gff3_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GFF3 attributes field (key=value;key=value format).

    Args:
        attr_string: The 9th column of a GFF3 line.

    Returns:
        Dictionary of attribute key-value pairs.
    """
    attributes: dict[str, str] = {}
    if not attr_string or attr_string == ".":
        return attributes

    for pair in attr_string.split(";"):
        pair = pair.strip()
        if not pair:
            continue
        if "=" in pair:
            key, value = pair.split("=", 1)
            # URL-decode common escapes
            value = value.replace("%3B", ";").replace("%3D", "=")
            value = value.replace("%26", "&").replace("%2C", ",")
            attributes[key] = value

    return attributes


def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF attributes field (key "value"; key "value"; format).

    Args:
        attr_string: The 9th column of a GTF line.

    Returns:
        Dictionary of attribute key-value pairs.
    """
    attributes: dict[str, str] = {}
    if not attr_string or attr_string == ".":
        return attributes

    # GTF format: key "value"; key "value";
    pattern = re.compile(r'(\w+)\s+"([^"]*)"')
    for match in pattern.finditer(attr_string):
        key, value = match.groups()
        attributes[key] = value

    return attributes


def validate_gff_line(
    line: str,
    line_num: int,
    is_gtf: bool = False
) -> tuple[bool, Optional[str]]:
    """
    Validate a GFF3/GTF line format.

    Args:
        line: The line to validate.
        line_num: Line number for error reporting.
        is_gtf: Whether to use GTF format rules.

    Returns:
        Tuple of (is_valid, error_message).
        If valid, error_message is None.
    """
    line = line.strip()

    # Skip empty lines and comments
    if not line or line.startswith("#"):
        return True, None

    fields = line.split("\t")

    # Check field count
    if len(fields) != 9:
        return False, f"Expected 9 tab-separated fields, got {len(fields)}"

    # Validate start/end are integers
    try:
        start = int(fields[3])
        end = int(fields[4])
    except ValueError:
        return False, f"Start ({fields[3]}) and end ({fields[4]}) must be integers"

    # Validate start <= end
    if start > end:
        return False, f"Start ({start}) cannot be greater than end ({end})"

    # Validate strand
    if fields[6] not in ("+", "-", "."):
        return False, f"Invalid strand '{fields[6]}', must be +, -, or ."

    # Validate phase for CDS
    if fields[2] == "CDS" and fields[7] not in ("0", "1", "2", "."):
        return False, f"Invalid phase '{fields[7]}' for CDS, must be 0, 1, 2, or ."

    return True, None


def parse_gff_line(
    line: str,
    line_num: int,
    is_gtf: bool = False
) -> Optional[GFFFeature]:
    """
    Parse a single GFF3/GTF line into a GFFFeature.

    Args:
        line: The line to parse.
        line_num: Line number for error reporting.
        is_gtf: Whether to use GTF format rules.

    Returns:
        GFFFeature object, or None for comments/empty lines.

    Raises:
        CompGeneError: If line format is invalid.
    """
    line = line.strip()

    # Skip empty lines and comments
    if not line or line.startswith("#"):
        return None

    # Validate line format
    is_valid, error_msg = validate_gff_line(line, line_num, is_gtf)
    if not is_valid:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"Invalid {'GTF' if is_gtf else 'GFF3'} format at line {line_num}",
            details=error_msg
        )

    fields = line.split("\t")

    # Parse score (can be ".")
    score: Optional[float] = None
    if fields[5] != ".":
        try:
            score = float(fields[5])
        except ValueError:
            pass

    # Parse phase (can be ".")
    phase: Optional[int] = None
    if fields[7] != ".":
        try:
            phase = int(fields[7])
        except ValueError:
            pass

    # Parse attributes based on format
    if is_gtf:
        attributes = parse_gtf_attributes(fields[8])
    else:
        attributes = parse_gff3_attributes(fields[8])

    return GFFFeature(
        seqid=fields[0],
        source=fields[1],
        type=fields[2],
        start=int(fields[3]),
        end=int(fields[4]),
        score=score,
        strand=fields[6],
        phase=phase,
        attributes=attributes,
        line_num=line_num
    )


def parse_gff3(path: Path) -> Iterator[GFFFeature]:
    """
    Stream-parse a GFF3 file, yielding features one at a time.

    Memory efficient - does not load entire file into memory.
    Supports gzip-compressed files (.gff3.gz).

    Args:
        path: Path to the GFF3 file.

    Yields:
        GFFFeature objects for each valid feature line.

    Raises:
        CompGeneError: If file not found or format is invalid.

    Example:
        >>> for feature in parse_gff3(Path("annotation.gff3")):
        ...     if feature.type == "gene":
        ...         print(feature.id)
    """
    with open_gff_file(path) as f:
        for line_num, line in enumerate(f, 1):
            feature = parse_gff_line(line, line_num, is_gtf=False)
            if feature is not None:
                yield feature


def parse_gtf(path: Path) -> Iterator[GFFFeature]:
    """
    Stream-parse a GTF file, yielding features one at a time.

    Memory efficient - does not load entire file into memory.
    Supports gzip-compressed files (.gtf.gz).

    Args:
        path: Path to the GTF file.

    Yields:
        GFFFeature objects for each valid feature line.

    Raises:
        CompGeneError: If file not found or format is invalid.

    Example:
        >>> for feature in parse_gtf(Path("annotation.gtf")):
        ...     if feature.type == "gene":
        ...         print(feature.attributes.get("gene_id"))
    """
    with open_gff_file(path) as f:
        for line_num, line in enumerate(f, 1):
            feature = parse_gff_line(line, line_num, is_gtf=True)
            if feature is not None:
                yield feature


# =============================================================================
# Gene Hierarchy Construction
# =============================================================================

# Feature types that represent transcripts
TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lncRNA", "rRNA", "tRNA", "snRNA", "snoRNA"}

# Feature types that represent coding sequences
CDS_TYPES = {"CDS"}

# Feature types that represent exons
EXON_TYPES = {"exon"}


def _get_parent_transcript_ids(feature: GFFFeature) -> list[str]:
    """
    Get parent transcript ID(s) for a CDS/exon feature.

    Handles both GFF3 (Parent attribute, possibly comma-separated) and
    GTF (transcript_id attribute) formats.

    Args:
        feature: A CDS or exon GFFFeature.

    Returns:
        List of parent transcript IDs.
    """
    # GTF format: use transcript_id attribute directly
    transcript_id = feature.attributes.get("transcript_id")
    if transcript_id:
        return [transcript_id]

    # GFF3 format: Parent attribute may be comma-separated for multi-parent
    parent = feature.attributes.get("Parent")
    if parent:
        return [p.strip() for p in parent.split(",") if p.strip()]

    return []


def build_gene_hierarchy(
    features: Iterator[GFFFeature]
) -> dict[str, GeneModel]:
    """
    Build gene → transcript → CDS hierarchy from parsed features.

    This function consumes the feature iterator and builds a complete
    gene hierarchy with transcripts and their exon/CDS features.

    Handles both GFF3 and GTF formats correctly:
    - GFF3: Uses Parent attribute (supports multi-parent features)
    - GTF: Uses transcript_id attribute for CDS/exon assignment

    Args:
        features: Iterator of GFFFeature objects (from parse_gff3/parse_gtf).

    Returns:
        Dictionary mapping gene_id to GeneModel objects.

    Example:
        >>> features = parse_gff3(Path("annotation.gff3"))
        >>> genes = build_gene_hierarchy(features)
        >>> for gene_id, gene in genes.items():
        ...     print(f"{gene_id}: {len(gene.transcripts)} transcripts")
    """
    genes: dict[str, GeneModel] = {}
    transcripts: dict[str, TranscriptModel] = {}

    # First pass: collect all features
    all_features: list[GFFFeature] = []
    for feature in features:
        all_features.append(feature)

        # Create gene entries
        if feature.type == "gene":
            gene_id = feature.id
            if gene_id:
                genes[gene_id] = GeneModel(gene_id=gene_id, feature=feature)

        # Create transcript entries
        elif feature.type in TRANSCRIPT_TYPES:
            transcript_id = feature.id
            # For GTF, transcript_id is in attributes; for GFF3, use ID
            if not transcript_id:
                transcript_id = feature.attributes.get("transcript_id")

            # Get parent gene: GFF3 uses Parent, GTF uses gene_id
            parent_gene = feature.attributes.get("Parent") or feature.attributes.get("gene_id")

            if transcript_id and parent_gene:
                transcripts[transcript_id] = TranscriptModel(
                    transcript_id=transcript_id,
                    gene_id=parent_gene,
                    feature=feature
                )

    # Second pass: assign exons and CDS to transcripts
    for feature in all_features:
        if feature.type in EXON_TYPES:
            parent_ids = _get_parent_transcript_ids(feature)
            for parent_id in parent_ids:
                if parent_id in transcripts:
                    transcripts[parent_id].exons.append(feature)
        elif feature.type in CDS_TYPES:
            parent_ids = _get_parent_transcript_ids(feature)
            for parent_id in parent_ids:
                if parent_id in transcripts:
                    transcripts[parent_id].cds_features.append(feature)

    # Third pass: assign transcripts to genes
    for transcript in transcripts.values():
        gene_id = transcript.gene_id
        if gene_id not in genes:
            # Create implicit gene if not found
            genes[gene_id] = GeneModel(gene_id=gene_id)
        genes[gene_id].transcripts[transcript.transcript_id] = transcript

    # Select representative transcripts
    for gene in genes.values():
        gene.representative_transcript = get_representative_transcript(gene)

    return genes


def get_representative_transcript(gene: GeneModel) -> Optional[str]:
    """
    Select the representative transcript for a gene.

    Selection criteria (in order):
    1. Longest CDS total length
    2. If no CDS, longest exon total length
    3. Tie-breaker: lexicographically smallest transcript ID

    Args:
        gene: GeneModel with transcripts.

    Returns:
        Transcript ID of the representative, or None if no transcripts.

    Example:
        >>> rep_id = get_representative_transcript(gene)
        >>> rep_transcript = gene.transcripts[rep_id]
    """
    if not gene.transcripts:
        return None

    # Sort by (effective_length descending, transcript_id ascending)
    sorted_transcripts = sorted(
        gene.transcripts.values(),
        key=lambda t: (-t.effective_length, t.transcript_id)
    )

    return sorted_transcripts[0].transcript_id


# =============================================================================
# Format Detection
# =============================================================================

def detect_gff_format(path: Path) -> str:
    """
    Detect whether a file is GFF3 or GTF format.

    Detection is based on:
    1. File extension (.gff3, .gff, .gtf)
    2. Header lines (##gff-version)
    3. Attribute format sampling

    Args:
        path: Path to the annotation file.

    Returns:
        "gff3" or "gtf" string.

    Raises:
        CompGeneError: If format cannot be determined.
    """
    # Check extension first
    name = path.name.lower()
    if name.endswith(".gff3") or name.endswith(".gff3.gz"):
        return "gff3"
    if name.endswith(".gtf") or name.endswith(".gtf.gz"):
        return "gtf"

    # Read first lines to determine format
    with open_gff_file(path) as f:
        for i, line in enumerate(f):
            if i > 100:  # Sample first 100 lines
                break

            line = line.strip()

            # Check for GFF3 version header
            if line.startswith("##gff-version"):
                return "gff3"

            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) == 9:
                attr = fields[8]
                # GTF uses key "value"; format
                if re.search(r'\w+\s+"[^"]*"', attr):
                    return "gtf"
                # GFF3 uses key=value; format
                if "=" in attr:
                    return "gff3"

    # Default to GFF3 if uncertain
    return "gff3"


def parse_annotation(path: Path) -> Iterator[GFFFeature]:
    """
    Parse an annotation file, auto-detecting GFF3 vs GTF format.

    This convenience function automatically detects the annotation format
    based on file extension and content, then delegates to the appropriate
    parser. Supports gzip-compressed files.

    Args:
        path: Path to the annotation file (.gff3, .gtf, or .gz variants).

    Yields:
        GFFFeature objects for each valid feature line.

    Raises:
        CompGeneError: If file not found or format is invalid.

    Example:
        >>> for feature in parse_annotation(Path("annotation.gff3")):
        ...     if feature.type == "gene":
        ...         print(f"Gene: {feature.id}")
        >>> # Also works with GTF
        >>> genes = build_gene_hierarchy(parse_annotation(Path("genes.gtf")))
    """
    fmt = detect_gff_format(path)
    if fmt == "gtf":
        yield from parse_gtf(path)
    else:
        yield from parse_gff3(path)


def validate_gff(path: Path, max_lines: int = 1000) -> tuple[bool, Optional[str]]:
    """
    Validate GFF3/GTF file format.

    Args:
        path: Path to the GFF/GTF file.
        max_lines: Maximum number of lines to check (default: 1000).

    Returns:
        Tuple of (is_valid, error_message).
        If valid, error_message is None.
    """
    try:
        fmt = detect_gff_format(path)
        is_gtf = (fmt == "gtf")

        with open_gff_file(path) as f:
            has_features = False

            for line_num, line in enumerate(f, 1):
                if line_num > max_lines:
                    break

                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                is_valid, error_msg = validate_gff_line(line, line_num, is_gtf)
                if not is_valid:
                    return False, f"Line {line_num}: {error_msg}"
                has_features = True

            if not has_features:
                return False, "No valid feature lines found"

            return True, None

    except CompGeneError as e:
        return False, str(e)


# =============================================================================
# Format Detection Utilities
# =============================================================================

def detect_format(path: Path) -> str:
    """
    Detect file format (gff3, gtf, or fasta).

    Args:
        path: Path to the file.

    Returns:
        Format string: "gff3", "gtf", or "fasta".

    Raises:
        CompGeneError: If format cannot be determined.
    """
    name = path.name.lower()

    # Check extension first
    if name.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz", ".fna", ".fna.gz")):
        return "fasta"
    if name.endswith((".gff3", ".gff3.gz", ".gff", ".gff.gz")):
        return "gff3"
    if name.endswith((".gtf", ".gtf.gz")):
        return "gtf"

    # Sample file content
    try:
        with open_gff_file(path) as f:
            for i, line in enumerate(f):
                if i > 20:
                    break

                line = line.strip()
                if not line:
                    continue

                # FASTA starts with >
                if line.startswith(">"):
                    return "fasta"

                # GFF3 version header
                if line.startswith("##gff-version"):
                    return "gff3"

                # Skip comments
                if line.startswith("#"):
                    continue

                # Check for tab-separated format
                fields = line.split("\t")
                if len(fields) == 9:
                    attr = fields[8]
                    if re.search(r'\w+\s+"[^"]*"', attr):
                        return "gtf"
                    if "=" in attr:
                        return "gff3"
    except Exception:
        pass

    raise CompGeneError(
        ErrorCode.E_INPUT_FORMAT,
        f"Cannot determine format of file: {path}",
        details="File does not appear to be GFF3, GTF, or FASTA format"
    )


def validate_file(path: Path, expected_format: Optional[str] = None) -> tuple[bool, Optional[str]]:
    """
    Validate a file's format.

    Args:
        path: Path to the file.
        expected_format: Expected format ("gff3", "gtf", "fasta"), or None to auto-detect.

    Returns:
        Tuple of (is_valid, error_message).
        If valid, error_message is None.
    """
    # Import here to avoid circular imports
    from workflow.lib.fasta import validate_fasta

    try:
        detected = detect_format(path)

        if expected_format and detected != expected_format:
            return False, f"Expected {expected_format} format, detected {detected}"

        if detected == "fasta":
            return validate_fasta(path)
        else:
            return validate_gff(path)

    except CompGeneError as e:
        return False, str(e)
