"""
CompGene FASTA Parser Module.

This module provides streaming parsers for FASTA sequence files,
with support for reading, writing, and computing sequence statistics.

Key features:
- Streaming/iterator-based parsing for memory efficiency (NFR2)
- gzip compression support
- Atomic file writing for checkpoint safety
- Sequence statistics computation

Source: Story 2.1 - GFF/FASTA 解析库
"""

import gzip
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, Iterator, Optional

from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.io import atomic_write


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class FastaRecord:
    """
    Represents a single sequence record from a FASTA file.

    Attributes:
        id: Sequence identifier (first word after ">").
        description: Full description line (everything after ">").
        sequence: The nucleotide or amino acid sequence.
    """

    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        """Get the sequence length in characters."""
        return len(self.sequence)

    def __len__(self) -> int:
        """Support len() on FastaRecord."""
        return self.length


@dataclass
class SequenceStats:
    """
    Statistics for a collection of sequences.

    Attributes:
        count: Number of sequences.
        total_length: Total length of all sequences.
        min_length: Minimum sequence length.
        max_length: Maximum sequence length.
        mean_length: Mean sequence length.
        n50: N50 value (length at which 50% of total bases are in longer sequences).
        gc_content: GC content as a fraction (0.0 to 1.0), or None for proteins.
    """

    count: int
    total_length: int
    min_length: int
    max_length: int
    mean_length: float
    n50: int
    gc_content: Optional[float] = None


# =============================================================================
# File Utilities
# =============================================================================

def open_fasta_file(path: Path, mode: str = "rt") -> IO[str]:
    """
    Open a FASTA file, automatically detecting gzip compression.

    Args:
        path: Path to the file.
        mode: File open mode (default: "rt" for text reading).

    Returns:
        File handle for reading/writing.

    Raises:
        CompGeneError: If file cannot be opened.
    """
    if "r" in mode and not path.exists():
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"FASTA file not found: {path}",
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
# FASTA Parsing
# =============================================================================

def parse_fasta(path: Path) -> Iterator[FastaRecord]:
    """
    Stream-parse a FASTA file, yielding records one at a time.

    Memory efficient - does not load entire file into memory.
    Supports gzip-compressed files (.fa.gz, .fasta.gz).

    Args:
        path: Path to the FASTA file.

    Yields:
        FastaRecord objects for each sequence.

    Raises:
        CompGeneError: If file not found or format is invalid.

    Example:
        >>> for record in parse_fasta(Path("sequences.fa")):
        ...     print(f"{record.id}: {record.length} bp")
    """
    current_id: Optional[str] = None
    current_description: Optional[str] = None
    current_sequence: list[str] = []
    line_num = 0

    with open_fasta_file(path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip("\n\r")

            if line.startswith(">"):
                # Yield previous record if exists
                if current_id is not None:
                    yield FastaRecord(
                        id=current_id,
                        description=current_description or current_id,
                        sequence="".join(current_sequence)
                    )

                # Parse header line
                header = line[1:]  # Remove ">"
                parts = header.split(None, 1)  # Split on first whitespace
                if not parts:
                    raise CompGeneError(
                        ErrorCode.E_INPUT_FORMAT,
                        f"Empty FASTA header at line {line_num}",
                        details="Header line must contain at least a sequence ID"
                    )
                current_id = parts[0]
                current_description = header
                current_sequence = []

            elif current_id is not None:
                # Sequence line - strip and add
                current_sequence.append(line.strip())

            elif line.strip():
                # Non-empty line before any header
                raise CompGeneError(
                    ErrorCode.E_INPUT_FORMAT,
                    f"Sequence data before header at line {line_num}",
                    details="FASTA file must start with a header line (>)"
                )

    # Yield final record
    if current_id is not None:
        yield FastaRecord(
            id=current_id,
            description=current_description or current_id,
            sequence="".join(current_sequence)
        )


def parse_fasta_to_dict(path: Path) -> dict[str, FastaRecord]:
    """
    Parse a FASTA file into a dictionary keyed by sequence ID.

    Note: Loads entire file into memory. For large files, use parse_fasta().

    Args:
        path: Path to the FASTA file.

    Returns:
        Dictionary mapping sequence ID to FastaRecord.

    Raises:
        CompGeneError: If file not found, format invalid, or duplicate IDs.
    """
    records: dict[str, FastaRecord] = {}
    for record in parse_fasta(path):
        if record.id in records:
            raise CompGeneError(
                ErrorCode.E_INPUT_FORMAT,
                f"Duplicate sequence ID: {record.id}",
                details="FASTA file contains duplicate sequence identifiers"
            )
        records[record.id] = record
    return records


# =============================================================================
# FASTA Writing
# =============================================================================

def format_fasta_record(record: FastaRecord, line_width: int = 60) -> str:
    """
    Format a FastaRecord as FASTA text.

    Args:
        record: The FastaRecord to format.
        line_width: Maximum line width for sequence (default: 60).

    Returns:
        Formatted FASTA string including header and sequence.
    """
    lines = [f">{record.description}"]

    # Wrap sequence at specified width
    seq = record.sequence
    for i in range(0, len(seq), line_width):
        lines.append(seq[i:i + line_width])

    return "\n".join(lines)


def write_fasta(
    records: Iterable[FastaRecord],
    path: Path,
    line_width: int = 60
) -> int:
    """
    Write FASTA records to a file using atomic write.

    Uses atomic write pattern for checkpoint safety - writes to a temporary
    file first, then renames to the target path.

    Args:
        records: Iterable of FastaRecord objects.
        path: Output file path.
        line_width: Maximum line width for sequence (default: 60).

    Returns:
        Number of records written.

    Raises:
        OSError: If write fails.

    Example:
        >>> records = [FastaRecord("seq1", "seq1 description", "ATCG")]
        >>> write_fasta(records, Path("output.fa"))
    """
    # Collect all content first
    content_lines: list[str] = []
    count = 0

    for record in records:
        content_lines.append(format_fasta_record(record, line_width))
        count += 1

    content = "\n".join(content_lines)
    if content:
        content += "\n"  # Trailing newline

    atomic_write(path, content)
    return count


def write_fasta_gzip(
    records: Iterable[FastaRecord],
    path: Path,
    line_width: int = 60
) -> int:
    """
    Write FASTA records to a gzip-compressed file using atomic write.

    Uses atomic write pattern for checkpoint safety - writes to a temporary
    file first, then renames to the target path.

    Args:
        records: Iterable of FastaRecord objects.
        path: Output file path (should end with .gz).
        line_width: Maximum line width for sequence (default: 60).

    Returns:
        Number of records written.

    Raises:
        OSError: If write or rename fails (after cleanup).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    count = 0

    try:
        with gzip.open(temp_path, "wt", encoding="utf-8") as f:
            for record in records:
                f.write(format_fasta_record(record, line_width))
                f.write("\n")
                count += 1

        # Atomic rename
        os.rename(temp_path, path)
    except Exception:
        # Clean up temporary file on failure
        if temp_path.exists():
            temp_path.unlink()
        raise

    return count


# =============================================================================
# Sequence Statistics
# =============================================================================

def get_sequence_length(record: FastaRecord) -> int:
    """
    Get the length of a sequence.

    Args:
        record: FastaRecord object.

    Returns:
        Sequence length in characters.
    """
    return len(record.sequence)


def compute_gc_content(sequence: str) -> Optional[float]:
    """
    Compute GC content of a nucleotide sequence.

    Args:
        sequence: Nucleotide sequence string.

    Returns:
        GC content as a fraction (0.0 to 1.0), or None if sequence is empty
        or appears to be protein.
    """
    if not sequence:
        return None

    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    atgc_count = (
        sequence.count("A") +
        sequence.count("T") +
        sequence.count("G") +
        sequence.count("C")
    )

    # If less than 80% ATGC, likely protein sequence
    if atgc_count < len(sequence) * 0.8:
        return None

    if atgc_count == 0:
        return None

    return gc_count / atgc_count


def compute_n50(lengths: list[int]) -> int:
    """
    Compute N50 value for a list of sequence lengths.

    N50 is the length at which 50% of the total sequence
    is contained in sequences of that length or longer.

    Args:
        lengths: List of sequence lengths.

    Returns:
        N50 value.
    """
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    half = total / 2
    cumsum = 0

    for length in sorted_lengths:
        cumsum += length
        if cumsum >= half:
            return length

    return sorted_lengths[-1]


def get_sequence_stats(records: Iterable[FastaRecord]) -> SequenceStats:
    """
    Compute statistics for a collection of sequences.

    Args:
        records: Iterable of FastaRecord objects.

    Returns:
        SequenceStats object with computed statistics.

    Example:
        >>> stats = get_sequence_stats(parse_fasta(Path("sequences.fa")))
        >>> print(f"N50: {stats.n50}, Count: {stats.count}")
    """
    lengths: list[int] = []
    total_gc_bases = 0
    total_atgc_bases = 0

    for record in records:
        lengths.append(record.length)

        # Compute GC for nucleotide sequences
        seq = record.sequence.upper()
        gc = seq.count("G") + seq.count("C")
        atgc = seq.count("A") + seq.count("T") + gc
        total_gc_bases += gc
        total_atgc_bases += atgc

    if not lengths:
        return SequenceStats(
            count=0,
            total_length=0,
            min_length=0,
            max_length=0,
            mean_length=0.0,
            n50=0,
            gc_content=None
        )

    total_length = sum(lengths)
    gc_content: Optional[float] = None
    if total_atgc_bases > total_length * 0.8:  # Likely nucleotide
        gc_content = total_gc_bases / total_atgc_bases if total_atgc_bases > 0 else None

    return SequenceStats(
        count=len(lengths),
        total_length=total_length,
        min_length=min(lengths),
        max_length=max(lengths),
        mean_length=total_length / len(lengths),
        n50=compute_n50(lengths),
        gc_content=gc_content
    )


# =============================================================================
# Validation
# =============================================================================

def validate_fasta(path: Path, max_lines: int = 1000) -> tuple[bool, Optional[str]]:
    """
    Validate FASTA file format.

    Args:
        path: Path to the FASTA file.
        max_lines: Maximum number of lines to check (default: 1000).

    Returns:
        Tuple of (is_valid, error_message).
        If valid, error_message is None.
    """
    try:
        with open_fasta_file(path) as f:
            has_header = False
            has_sequence = False
            current_has_sequence = False
            has_content = False

            for line_num, line in enumerate(f, 1):
                if line_num > max_lines:
                    break

                line = line.strip()
                if not line:
                    continue

                has_content = True

                if line.startswith(">"):
                    if has_header and not current_has_sequence:
                        return False, f"Empty sequence at line {line_num - 1}"
                    has_header = True
                    current_has_sequence = False
                elif has_header:
                    has_sequence = True
                    current_has_sequence = True
                else:
                    return False, f"Sequence data before header at line {line_num}"

            # Check for empty or whitespace-only files
            if not has_content:
                return False, "Empty file or file contains only whitespace"

            # Check for file with no headers
            if not has_header:
                return False, "No FASTA headers found in file"

            if has_header and not has_sequence:
                return False, "No sequence data found after header"

            return True, None

    except CompGeneError as e:
        return False, str(e)
