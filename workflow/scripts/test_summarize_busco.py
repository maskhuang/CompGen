"""
Unit tests for workflow/scripts/summarize_busco.py

Tests for BUSCO summary parsing and aggregation.
"""

from pathlib import Path
import csv
import pytest

# Import the module under test
from workflow.scripts.summarize_busco import (
    parse_summary,
    BUSCO_SUMMARY_PATTERN,
    LINEAGE_PATTERN,
    VERSION_PATTERN,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def sample_summary_content() -> str:
    """Sample BUSCO short_summary.txt content."""
    return """# BUSCO version is: 5.4.7
# The lineage dataset is: primates_odb10 (Creation date: 2024-01-08, number of BUSCOs: 13780)
# Summarized benchmarking in BUSCO notation for file proteins.fa
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:13780
        13118   Complete BUSCOs (C)
        12828   Complete and single-copy BUSCOs (S)
        290     Complete and duplicated BUSCOs (D)
        317     Fragmented BUSCOs (F)
        345     Missing BUSCOs (M)
        13780   Total BUSCO groups searched
"""


@pytest.fixture
def busco6_summary_content() -> str:
    """Sample BUSCO 6.x short_summary.txt content."""
    return """# BUSCO version is: 6.0.0
# The lineage dataset is: eukaryota_odb10 (Creation date: 2024-06-01, number of BUSCOs: 255)
# Summarized benchmarking in BUSCO notation for file proteins.fa
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:98.4%[S:96.1%,D:2.3%],F:1.2%,M:0.4%,n:255
        251     Complete BUSCOs (C)
        245     Complete and single-copy BUSCOs (S)
        6       Complete and duplicated BUSCOs (D)
        3       Fragmented BUSCOs (F)
        1       Missing BUSCOs (M)
        255     Total BUSCO groups searched
"""


# =============================================================================
# Test parse_summary
# =============================================================================

class TestParseSummary:
    """Tests for parse_summary function."""

    def test_parse_summary_busco5(self, tmp_path: Path, sample_summary_content: str) -> None:
        """[P1] Given BUSCO 5.x summary file, when parsed, then correct values extracted."""
        # Given
        # Create directory structure: results/qc/{species}/busco/short_summary.txt
        species_dir = tmp_path / "results" / "qc" / "mmur" / "busco"
        species_dir.mkdir(parents=True)
        summary_file = species_dir / "short_summary.txt"
        summary_file.write_text(sample_summary_content)

        # When
        result = parse_summary(summary_file)

        # Then
        assert result["species"] == "mmur"
        assert result["complete_pct"] == 95.2
        assert result["single_copy_pct"] == 93.1
        assert result["duplicated_pct"] == 2.1
        assert result["fragmented_pct"] == 2.3
        assert result["missing_pct"] == 2.5
        assert result["total"] == 13780
        assert result["lineage"] == "primates_odb10"
        assert result["busco_version"] == "5.4.7"

    def test_parse_summary_busco6(self, tmp_path: Path, busco6_summary_content: str) -> None:
        """[P1] Given BUSCO 6.x summary file, when parsed, then correct values extracted."""
        # Given
        species_dir = tmp_path / "results" / "qc" / "hsap" / "busco"
        species_dir.mkdir(parents=True)
        summary_file = species_dir / "short_summary.txt"
        summary_file.write_text(busco6_summary_content)

        # When
        result = parse_summary(summary_file)

        # Then
        assert result["species"] == "hsap"
        assert result["complete_pct"] == 98.4
        assert result["single_copy_pct"] == 96.1
        assert result["duplicated_pct"] == 2.3
        assert result["fragmented_pct"] == 1.2
        assert result["missing_pct"] == 0.4
        assert result["total"] == 255
        assert result["lineage"] == "eukaryota_odb10"
        assert result["busco_version"] == "6.0.0"

    def test_parse_summary_invalid_content_raises(self, tmp_path: Path) -> None:
        """[P1] Given invalid content, when parsed, then ValueError raised."""
        # Given
        summary_file = tmp_path / "invalid.txt"
        summary_file.write_text("Invalid content without BUSCO results")

        # When/Then
        with pytest.raises(ValueError) as exc_info:
            parse_summary(summary_file)
        assert "parse" in str(exc_info.value).lower() or "busco" in str(exc_info.value).lower()

    def test_parse_summary_extracts_species_from_path(self, tmp_path: Path, sample_summary_content: str) -> None:
        """[P2] Given summary file, when parsed, then species extracted from path."""
        # Given - species name is in grandparent directory
        species_dir = tmp_path / "macaque_01" / "busco"
        species_dir.mkdir(parents=True)
        summary_file = species_dir / "short_summary.txt"
        summary_file.write_text(sample_summary_content)

        # When
        result = parse_summary(summary_file)

        # Then
        assert result["species"] == "macaque_01"


# =============================================================================
# Test regex patterns
# =============================================================================

class TestRegexPatterns:
    """Tests for regex patterns used in parsing."""

    def test_busco_summary_pattern_matches_notation(self) -> None:
        """[P1] Given BUSCO notation, when pattern applied, then matches."""
        # Given
        notation = "C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:13780"

        # When
        match = BUSCO_SUMMARY_PATTERN.search(notation)

        # Then
        assert match is not None
        assert match.group(1) == "95.2"
        assert match.group(2) == "93.1"
        assert match.group(3) == "2.1"
        assert match.group(4) == "2.3"
        assert match.group(5) == "2.5"
        assert match.group(6) == "13780"

    def test_busco_summary_pattern_handles_integers(self) -> None:
        """[P2] Given integer percentages, when pattern applied, then matches."""
        # Given - some BUSCO outputs use integers (e.g., 100%)
        notation = "C:100%[S:98%,D:2%],F:0%,M:0%,n:255"

        # When
        match = BUSCO_SUMMARY_PATTERN.search(notation)

        # Then
        assert match is not None
        assert match.group(1) == "100"
        assert match.group(5) == "0"

    def test_lineage_pattern_extracts_lineage(self) -> None:
        """[P1] Given lineage line, when pattern applied, then lineage extracted."""
        # Given
        line = "# The lineage dataset is: primates_odb10 (Creation date: 2024-01-08)"

        # When
        match = LINEAGE_PATTERN.search(line)

        # Then
        assert match is not None
        assert match.group(1) == "primates_odb10"

    def test_version_pattern_extracts_version(self) -> None:
        """[P1] Given version line, when pattern applied, then version extracted."""
        # Given
        line = "# BUSCO version is: 5.4.7"

        # When
        match = VERSION_PATTERN.search(line)

        # Then
        assert match is not None
        assert match.group(1) == "5.4.7"
