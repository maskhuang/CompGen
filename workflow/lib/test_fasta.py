"""
Unit tests for the FASTA parser module.

Tests cover:
- FASTA parsing and record extraction
- Sequence statistics computation
- FASTA writing with atomic operations
- gzip compression support
- Error handling for malformed input
- Memory efficiency (streaming)

Source: Story 2.1 - GFF/FASTA 解析库
"""

import gzip
import tempfile
from pathlib import Path

import pytest

from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.fasta import (
    FastaRecord,
    SequenceStats,
    parse_fasta,
    parse_fasta_to_dict,
    write_fasta,
    write_fasta_gzip,
    format_fasta_record,
    get_sequence_length,
    get_sequence_stats,
    compute_gc_content,
    compute_n50,
    validate_fasta,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def fasta_content() -> str:
    """Sample FASTA content for testing."""
    return """>seq1 First sequence
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCG
>seq2 Second sequence with longer header description
GCTAGCTAGCTAGCTAGCTAGCTAGCTA
>seq3 Protein-like sequence
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
"""


@pytest.fixture
def fasta_file(fasta_content: str, tmp_path: Path) -> Path:
    """Create a temporary FASTA file."""
    path = tmp_path / "test.fa"
    path.write_text(fasta_content)
    return path


@pytest.fixture
def fasta_gzip_file(fasta_content: str, tmp_path: Path) -> Path:
    """Create a temporary gzip-compressed FASTA file."""
    path = tmp_path / "test.fa.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(fasta_content)
    return path


# =============================================================================
# Test FASTA Parsing
# =============================================================================

class TestFastaParsing:
    """Tests for FASTA file parsing."""

    def test_parse_fasta_file(self, fasta_file: Path):
        """Test parsing a FASTA file."""
        records = list(parse_fasta(fasta_file))

        assert len(records) == 3

        # Check first record
        assert records[0].id == "seq1"
        assert records[0].description == "seq1 First sequence"
        assert len(records[0].sequence) == 56  # 36 + 20

        # Check second record
        assert records[1].id == "seq2"
        assert "longer header" in records[1].description

    def test_parse_fasta_gzip(self, fasta_gzip_file: Path):
        """Test parsing a gzip-compressed FASTA file."""
        records = list(parse_fasta(fasta_gzip_file))
        assert len(records) == 3

    def test_parse_fasta_missing_file(self, tmp_path: Path):
        """Test error handling for missing file."""
        with pytest.raises(CompGeneError) as exc_info:
            list(parse_fasta(tmp_path / "nonexistent.fa"))
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_parse_fasta_empty_header(self, tmp_path: Path):
        """Test error handling for empty header."""
        bad_file = tmp_path / "bad.fa"
        bad_file.write_text(">\nATCG")

        with pytest.raises(CompGeneError) as exc_info:
            list(parse_fasta(bad_file))
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT

    def test_parse_fasta_no_header(self, tmp_path: Path):
        """Test error handling for sequence before header."""
        bad_file = tmp_path / "bad.fa"
        bad_file.write_text("ATCGATCG\n>seq1\nATCG")

        with pytest.raises(CompGeneError) as exc_info:
            list(parse_fasta(bad_file))
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT

    def test_parse_fasta_to_dict(self, fasta_file: Path):
        """Test parsing FASTA to dictionary."""
        records = parse_fasta_to_dict(fasta_file)

        assert "seq1" in records
        assert "seq2" in records
        assert "seq3" in records

    def test_parse_fasta_duplicate_ids(self, tmp_path: Path):
        """Test error handling for duplicate sequence IDs."""
        bad_file = tmp_path / "dup.fa"
        bad_file.write_text(">seq1\nATCG\n>seq1\nGCTA")

        with pytest.raises(CompGeneError) as exc_info:
            parse_fasta_to_dict(bad_file)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


# =============================================================================
# Test FastaRecord
# =============================================================================

class TestFastaRecord:
    """Tests for FastaRecord dataclass."""

    def test_record_length(self):
        """Test length property."""
        record = FastaRecord(id="seq1", description="seq1", sequence="ATCGATCG")
        assert record.length == 8
        assert len(record) == 8

    def test_record_empty_sequence(self):
        """Test record with empty sequence."""
        record = FastaRecord(id="empty", description="empty", sequence="")
        assert record.length == 0


# =============================================================================
# Test FASTA Writing
# =============================================================================

class TestFastaWriting:
    """Tests for FASTA file writing."""

    def test_write_fasta(self, tmp_path: Path):
        """Test writing FASTA file."""
        records = [
            FastaRecord("seq1", "seq1 description", "ATCGATCG"),
            FastaRecord("seq2", "seq2 another", "GCTAGCTA"),
        ]

        output_path = tmp_path / "output.fa"
        count = write_fasta(records, output_path)

        assert count == 2
        assert output_path.exists()

        # Read back and verify
        content = output_path.read_text()
        assert ">seq1 description" in content
        assert "ATCGATCG" in content

    def test_write_fasta_line_wrapping(self, tmp_path: Path):
        """Test line wrapping in FASTA output."""
        long_seq = "A" * 100
        records = [FastaRecord("seq1", "seq1", long_seq)]

        output_path = tmp_path / "output.fa"
        write_fasta(records, output_path, line_width=60)

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Should have header + 2 sequence lines (60 + 40)
        assert len(lines) == 3
        assert len(lines[1]) == 60
        assert len(lines[2]) == 40

    def test_write_fasta_gzip(self, tmp_path: Path):
        """Test writing gzip-compressed FASTA."""
        records = [FastaRecord("seq1", "seq1", "ATCGATCG")]

        output_path = tmp_path / "output.fa.gz"
        count = write_fasta_gzip(records, output_path)

        assert count == 1
        assert output_path.exists()

        # Read back and verify
        with gzip.open(output_path, "rt") as f:
            content = f.read()
        assert ">seq1" in content
        assert "ATCGATCG" in content

    def test_format_fasta_record(self):
        """Test formatting a single record."""
        record = FastaRecord("seq1", "seq1 description", "ATCGATCGATCG")
        formatted = format_fasta_record(record, line_width=8)

        lines = formatted.split("\n")
        assert lines[0] == ">seq1 description"
        assert lines[1] == "ATCGATCG"
        assert lines[2] == "ATCG"


# =============================================================================
# Test Sequence Statistics
# =============================================================================

class TestSequenceStats:
    """Tests for sequence statistics computation."""

    def test_get_sequence_length(self):
        """Test get_sequence_length function."""
        record = FastaRecord("seq1", "seq1", "ATCGATCG")
        assert get_sequence_length(record) == 8

    def test_compute_gc_content(self):
        """Test GC content computation."""
        # 50% GC
        gc = compute_gc_content("ATCGATCG")
        assert gc == pytest.approx(0.5, rel=0.01)

        # 100% GC
        gc = compute_gc_content("GGCC")
        assert gc == pytest.approx(1.0)

        # 0% GC
        gc = compute_gc_content("AATT")
        assert gc == pytest.approx(0.0)

        # Empty sequence
        gc = compute_gc_content("")
        assert gc is None

        # Protein sequence (should return None)
        gc = compute_gc_content("MVLSPADKTNVK")
        assert gc is None

    def test_compute_n50(self):
        """Test N50 computation."""
        # Simple case: [100, 50, 30, 20]
        # Total = 200, half = 100
        # Cumsum: 100 >= 100 at first element
        n50 = compute_n50([100, 50, 30, 20])
        assert n50 == 100

        # Another case: [60, 50, 40, 30, 20]
        # Total = 200, half = 100
        # Cumsum: 60, 110 >= 100 at second element
        n50 = compute_n50([60, 50, 40, 30, 20])
        assert n50 == 50

        # Empty list
        n50 = compute_n50([])
        assert n50 == 0

    def test_get_sequence_stats(self, fasta_file: Path):
        """Test full statistics computation."""
        records = parse_fasta(fasta_file)
        stats = get_sequence_stats(records)

        assert stats.count == 3
        assert stats.total_length > 0
        assert stats.min_length > 0
        assert stats.max_length > 0
        assert stats.mean_length > 0
        assert stats.n50 > 0

    def test_get_sequence_stats_empty(self):
        """Test statistics for empty input."""
        stats = get_sequence_stats([])

        assert stats.count == 0
        assert stats.total_length == 0
        assert stats.min_length == 0
        assert stats.max_length == 0
        assert stats.mean_length == 0.0
        assert stats.n50 == 0

    def test_gc_content_in_stats(self):
        """Test GC content in statistics."""
        records = [
            FastaRecord("s1", "s1", "ATCGATCG"),
            FastaRecord("s2", "s2", "GGCCGGCC"),
        ]
        stats = get_sequence_stats(records)

        # Both are nucleotide, average GC should be computed
        assert stats.gc_content is not None
        assert 0.0 <= stats.gc_content <= 1.0


# =============================================================================
# Test Validation
# =============================================================================

class TestValidation:
    """Tests for FASTA file validation."""

    def test_validate_fasta_valid(self, fasta_file: Path):
        """Test validation of valid FASTA file."""
        is_valid, error = validate_fasta(fasta_file)
        assert is_valid
        assert error is None

    def test_validate_fasta_no_header(self, tmp_path: Path):
        """Test validation of file without header."""
        bad_file = tmp_path / "bad.fa"
        bad_file.write_text("ATCGATCG")

        is_valid, error = validate_fasta(bad_file)
        assert not is_valid
        assert error is not None

    def test_validate_fasta_empty_sequence(self, tmp_path: Path):
        """Test validation of file with empty sequence."""
        bad_file = tmp_path / "bad.fa"
        bad_file.write_text(">seq1\n>seq2\nATCG")

        is_valid, error = validate_fasta(bad_file)
        assert not is_valid
        assert "Empty sequence" in error

    def test_validate_fasta_missing_file(self, tmp_path: Path):
        """Test validation of missing file."""
        is_valid, error = validate_fasta(tmp_path / "nonexistent.fa")
        assert not is_valid
        assert error is not None


# =============================================================================
# Test Streaming (Memory Efficiency)
# =============================================================================

class TestStreaming:
    """Tests for memory-efficient streaming."""

    def test_parse_fasta_is_generator(self, fasta_file: Path):
        """Test that parse_fasta returns a generator."""
        result = parse_fasta(fasta_file)
        assert hasattr(result, "__iter__")
        assert hasattr(result, "__next__")

    def test_parse_fasta_lazy_evaluation(self, tmp_path: Path):
        """Test that parsing is lazy (doesn't load all at once)."""
        # Create a file with multiple records
        content = "\n".join([f">seq{i}\nATCG" for i in range(100)])
        path = tmp_path / "many.fa"
        path.write_text(content)

        # Create generator
        gen = parse_fasta(path)

        # Get just first record
        first = next(gen)
        assert first.id == "seq0"

        # Generator should still have more
        second = next(gen)
        assert second.id == "seq1"


# =============================================================================
# Test Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_single_line_sequence(self, tmp_path: Path):
        """Test parsing single-line sequence."""
        path = tmp_path / "single.fa"
        path.write_text(">seq1\nATCG")

        records = list(parse_fasta(path))
        assert len(records) == 1
        assert records[0].sequence == "ATCG"

    def test_multiline_sequence(self, tmp_path: Path):
        """Test parsing multi-line sequence."""
        path = tmp_path / "multi.fa"
        path.write_text(">seq1\nATCG\nGCTA\nTTTT")

        records = list(parse_fasta(path))
        assert len(records) == 1
        assert records[0].sequence == "ATCGGCTATTTT"

    def test_empty_file(self, tmp_path: Path):
        """Test parsing empty file."""
        path = tmp_path / "empty.fa"
        path.write_text("")

        records = list(parse_fasta(path))
        assert len(records) == 0

    def test_only_whitespace(self, tmp_path: Path):
        """Test file with only whitespace returns empty list."""
        path = tmp_path / "whitespace.fa"
        path.write_text("\n\n   \n\n")

        records = list(parse_fasta(path))
        assert len(records) == 0

    def test_only_whitespace_validation_fails(self, tmp_path: Path):
        """Test that validation fails for whitespace-only file."""
        path = tmp_path / "whitespace.fa"
        path.write_text("\n\n   \n\n")

        is_valid, error = validate_fasta(path)
        assert not is_valid
        assert "Empty file" in error or "whitespace" in error

    def test_windows_line_endings(self, tmp_path: Path):
        """Test handling Windows line endings."""
        path = tmp_path / "windows.fa"
        path.write_text(">seq1\r\nATCG\r\nGCTA\r\n")

        records = list(parse_fasta(path))
        assert len(records) == 1
        assert records[0].sequence == "ATCGGCTA"


# =============================================================================
# Test Actual Fixture Files
# =============================================================================

class TestFixtureFiles:
    """Tests that verify the actual fixture files work correctly."""

    @pytest.fixture
    def fixtures_dir(self) -> Path:
        """Get the path to the fixtures directory."""
        return Path(__file__).parent.parent.parent / "tests" / "fixtures"

    def test_mini_genome_fixture(self, fixtures_dir: Path):
        """Test parsing the actual mini_genome.fa fixture file."""
        path = fixtures_dir / "mini_genome.fa"
        assert path.exists(), f"Fixture file not found: {path}"

        records = list(parse_fasta(path))

        # Should have 3 sequences: chr1, chr2, protein1
        assert len(records) == 3

        # Verify sequence IDs
        ids = [r.id for r in records]
        assert "chr1" in ids
        assert "chr2" in ids
        assert "protein1" in ids

        # Verify chr1 has nucleotide sequence (multi-line)
        chr1 = next(r for r in records if r.id == "chr1")
        assert chr1.length > 0
        assert all(c in "ATCGN" for c in chr1.sequence.upper())

        # Verify protein1 has protein-like sequence
        protein1 = next(r for r in records if r.id == "protein1")
        assert protein1.length > 0

    def test_mini_genome_fixture_validation(self, fixtures_dir: Path):
        """Test that mini_genome.fa passes validation."""
        path = fixtures_dir / "mini_genome.fa"
        is_valid, error = validate_fasta(path)
        assert is_valid, f"Fixture validation failed: {error}"

    def test_mini_genome_fixture_stats(self, fixtures_dir: Path):
        """Test computing statistics on mini_genome.fa fixture."""
        path = fixtures_dir / "mini_genome.fa"
        records = parse_fasta(path)
        stats = get_sequence_stats(records)

        assert stats.count == 3
        assert stats.total_length > 0
        assert stats.min_length > 0
        assert stats.max_length >= stats.min_length
        assert stats.n50 > 0
