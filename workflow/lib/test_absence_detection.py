"""
Unit tests for the absence detection module.

Tests cover:
- Gene classification logic (present/uncertain/missing)
- Coverage/identity extraction from GFF attributes
- Threshold parameterization
- Edge cases and boundary conditions
- Full GFF classification workflow

Source: Story 5.2 - 映射质量分析与缺失检测
"""

from pathlib import Path

import pytest

from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def liftoff_gff_content() -> str:
    """Sample Liftoff output GFF3 with coverage attributes."""
    return """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene001;Name=BRCA1;coverage=0.95;sequence_ID=NM_007294
chr1\tliftoff\tmRNA\t1000\t2000\t.\t+\t.\tID=mrna001;Parent=gene001;coverage=0.95
chr1\tliftoff\tgene\t3000\t4000\t.\t+\t.\tID=gene002;Name=TP53;coverage=0.42;sequence_ID=NM_000546
chr1\tliftoff\tmRNA\t3000\t4000\t.\t+\t.\tID=mrna002;Parent=gene002;coverage=0.42
chr1\tliftoff\tgene\t5000\t6000\t.\t+\t.\tID=gene003;Name=FOXP2;coverage=0.10;sequence_ID=NM_014491
chr1\tliftoff\tmRNA\t5000\t6000\t.\t+\t.\tID=mrna003;Parent=gene003;coverage=0.10
"""


@pytest.fixture
def liftoff_gff_file(liftoff_gff_content: str, tmp_path: Path) -> Path:
    """Create a temporary Liftoff GFF3 file."""
    path = tmp_path / "lifted_annotation.gff3"
    path.write_text(liftoff_gff_content)
    return path


@pytest.fixture
def unmapped_content() -> str:
    """Sample unmapped features file."""
    return """gene004
gene005
"""


@pytest.fixture
def unmapped_file(unmapped_content: str, tmp_path: Path) -> Path:
    """Create a temporary unmapped features file."""
    path = tmp_path / "unmapped_features.txt"
    path.write_text(unmapped_content)
    return path


@pytest.fixture
def empty_unmapped_file(tmp_path: Path) -> Path:
    """Create an empty unmapped features file."""
    path = tmp_path / "unmapped_features.txt"
    path.write_text("")
    return path


# =============================================================================
# Test Gene Classification Logic
# =============================================================================

class TestClassifyGene:
    """Tests for the classify_gene function."""

    def test_present_high_coverage_identity(self):
        """Test gene classified as present with high coverage and identity."""
        from workflow.lib.absence_detection import classify_gene

        result = classify_gene(coverage=0.95, identity=0.92)
        assert result == "present"

    def test_present_at_threshold(self):
        """Test gene classified as present at exactly the threshold."""
        from workflow.lib.absence_detection import classify_gene

        result = classify_gene(coverage=0.50, identity=0.50)
        assert result == "present"

    def test_uncertain_low_coverage(self):
        """Test gene classified as uncertain with low coverage."""
        from workflow.lib.absence_detection import classify_gene

        result = classify_gene(coverage=0.42, identity=0.80)
        assert result == "uncertain"

    def test_uncertain_low_identity(self):
        """Test gene classified as uncertain with low identity."""
        from workflow.lib.absence_detection import classify_gene

        result = classify_gene(coverage=0.80, identity=0.42)
        assert result == "uncertain"

    def test_uncertain_both_low(self):
        """Test gene classified as uncertain with both metrics low but non-zero."""
        from workflow.lib.absence_detection import classify_gene

        result = classify_gene(coverage=0.30, identity=0.30)
        assert result == "uncertain"

    def test_missing_zero_coverage(self):
        """Test gene classified as missing with zero coverage."""
        from workflow.lib.absence_detection import classify_gene

        result = classify_gene(coverage=0.0, identity=0.0)
        assert result == "missing"

    def test_custom_thresholds_stricter(self):
        """Test classification with stricter custom thresholds."""
        from workflow.lib.absence_detection import classify_gene

        # With stricter thresholds (0.8), 0.70 coverage should be uncertain
        result = classify_gene(
            coverage=0.70, identity=0.70,
            min_coverage=0.80, min_identity=0.80
        )
        assert result == "uncertain"

    def test_custom_thresholds_looser(self):
        """Test classification with looser custom thresholds."""
        from workflow.lib.absence_detection import classify_gene

        # With looser thresholds (0.3), 0.35 should be present
        result = classify_gene(
            coverage=0.35, identity=0.35,
            min_coverage=0.30, min_identity=0.30
        )
        assert result == "present"

    def test_identity_defaults_to_coverage(self):
        """Test that identity defaults to coverage if not provided."""
        from workflow.lib.absence_detection import classify_gene

        # When identity is None, use coverage as identity
        result = classify_gene(coverage=0.60, identity=None)
        assert result == "present"


# =============================================================================
# Test Coverage/Identity Extraction
# =============================================================================

class TestExtractCoverageIdentity:
    """Tests for extracting coverage and identity from GFF attributes."""

    def test_extract_coverage_from_attributes(self):
        """Test extracting coverage from GFF attribute dictionary."""
        from workflow.lib.absence_detection import extract_coverage_identity

        attributes = {"ID": "gene001", "coverage": "0.95", "sequence_ID": "NM_001"}
        coverage, identity = extract_coverage_identity(attributes)

        assert coverage == 0.95
        # Identity defaults to coverage when not present
        assert identity == 0.95

    def test_extract_both_coverage_and_identity(self):
        """Test extracting both coverage and identity when both present."""
        from workflow.lib.absence_detection import extract_coverage_identity

        attributes = {"coverage": "0.90", "identity": "0.85"}
        coverage, identity = extract_coverage_identity(attributes)

        assert coverage == 0.90
        assert identity == 0.85

    def test_missing_coverage_returns_zero(self):
        """Test that missing coverage returns 0.0."""
        from workflow.lib.absence_detection import extract_coverage_identity

        attributes = {"ID": "gene001", "Name": "BRCA1"}
        coverage, identity = extract_coverage_identity(attributes)

        assert coverage == 0.0
        assert identity == 0.0

    def test_invalid_coverage_value(self):
        """Test handling of invalid coverage value."""
        from workflow.lib.absence_detection import extract_coverage_identity

        attributes = {"coverage": "invalid"}
        coverage, identity = extract_coverage_identity(attributes)

        assert coverage == 0.0
        assert identity == 0.0

    def test_sequence_identity_attribute(self):
        """Test extracting sequence_identity attribute (alternative name)."""
        from workflow.lib.absence_detection import extract_coverage_identity

        attributes = {"coverage": "0.90", "sequence_identity": "0.88"}
        coverage, identity = extract_coverage_identity(attributes)

        assert coverage == 0.90
        assert identity == 0.88


# =============================================================================
# Test GFF Classification
# =============================================================================

class TestClassifyLiftoffGenes:
    """Tests for classifying genes from Liftoff GFF output."""

    def test_classify_genes_from_gff(self, liftoff_gff_file: Path, empty_unmapped_file: Path):
        """Test classifying genes from a Liftoff GFF file."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=empty_unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        # Should have 3 genes from GFF
        assert len(results) == 3

        # Check classifications
        gene001 = next(r for r in results if r["gene_id"] == "gene001")
        assert gene001["status"] == "present"
        assert gene001["coverage"] == 0.95

        gene002 = next(r for r in results if r["gene_id"] == "gene002")
        assert gene002["status"] == "uncertain"
        assert gene002["coverage"] == 0.42

        gene003 = next(r for r in results if r["gene_id"] == "gene003")
        assert gene003["status"] == "uncertain"  # coverage=0.10, below threshold
        assert gene003["coverage"] == 0.10

    def test_classify_with_unmapped_genes(
        self, liftoff_gff_file: Path, unmapped_file: Path
    ):
        """Test that unmapped genes are classified as missing."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        # Should have 5 genes total (3 from GFF + 2 unmapped)
        assert len(results) == 5

        # Check unmapped genes are missing
        gene004 = next(r for r in results if r["gene_id"] == "gene004")
        assert gene004["status"] == "missing"
        assert gene004["coverage"] == 0.0

        gene005 = next(r for r in results if r["gene_id"] == "gene005")
        assert gene005["status"] == "missing"

    def test_classify_with_custom_thresholds(
        self, liftoff_gff_file: Path, empty_unmapped_file: Path
    ):
        """Test classification with custom thresholds."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        # With very low threshold (0.1), gene003 should be present
        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=empty_unmapped_file,
            reference_species="human",
            target_species="mouse_lemur",
            min_coverage=0.10,
            min_identity=0.10
        )

        gene003 = next(r for r in results if r["gene_id"] == "gene003")
        assert gene003["status"] == "present"

    def test_result_contains_required_fields(
        self, liftoff_gff_file: Path, empty_unmapped_file: Path
    ):
        """Test that results contain all required fields."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=empty_unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        required_fields = [
            "gene_id", "gene_name", "reference_species", "target_species",
            "status", "coverage", "identity", "source_gene"
        ]

        for result in results:
            for field in required_fields:
                assert field in result, f"Missing field: {field}"

    def test_missing_gff_file_raises_error(self, tmp_path: Path, empty_unmapped_file: Path):
        """Test that missing GFF file raises appropriate error."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        with pytest.raises(CompGeneError) as exc_info:
            classify_liftoff_genes(
                lifted_gff=tmp_path / "nonexistent.gff3",
                unmapped_file=empty_unmapped_file,
                reference_species="human",
                target_species="mouse_lemur"
            )

        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_missing_unmapped_file_logs_warning(self, liftoff_gff_file: Path, tmp_path: Path, caplog):
        """Test that missing unmapped file logs warning but doesn't fail."""
        import logging
        from workflow.lib.absence_detection import classify_liftoff_genes

        # Use a non-existent unmapped file path
        nonexistent_unmapped = tmp_path / "does_not_exist.txt"

        with caplog.at_level(logging.WARNING):
            results = classify_liftoff_genes(
                lifted_gff=liftoff_gff_file,
                unmapped_file=nonexistent_unmapped,
                reference_species="human",
                target_species="mouse_lemur"
            )

        # Should still return results from GFF (3 genes)
        assert len(results) == 3

        # Should have logged a warning
        assert "Unmapped features file not found" in caplog.text


# =============================================================================
# Test Statistics Calculation
# =============================================================================

class TestCalculateClassificationStats:
    """Tests for calculating classification statistics."""

    def test_calculate_stats(self, liftoff_gff_file: Path, unmapped_file: Path):
        """Test calculating classification statistics."""
        from workflow.lib.absence_detection import (
            classify_liftoff_genes,
            calculate_classification_stats
        )

        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        stats = calculate_classification_stats(results)

        assert stats["total_genes"] == 5
        assert stats["present"] == 1  # Only gene001 has coverage >= 0.5
        assert stats["uncertain"] == 2  # gene002 and gene003 have 0 < coverage < 0.5
        assert stats["missing"] == 2  # gene004 and gene005 from unmapped
        assert stats["present_rate"] == pytest.approx(0.2)  # 1/5
        assert stats["missing_rate"] == pytest.approx(0.4)  # 2/5

    def test_stats_with_empty_results(self):
        """Test statistics calculation with empty results."""
        from workflow.lib.absence_detection import calculate_classification_stats

        stats = calculate_classification_stats([])

        assert stats["total_genes"] == 0
        assert stats["present"] == 0
        assert stats["uncertain"] == 0
        assert stats["missing"] == 0
        assert stats["present_rate"] == 0.0
        assert stats["missing_rate"] == 0.0


# =============================================================================
# Test Missing Genes Filter
# =============================================================================

class TestFilterMissingGenes:
    """Tests for filtering missing genes from classification results."""

    def test_filter_missing_genes(self, liftoff_gff_file: Path, unmapped_file: Path):
        """Test filtering only missing genes."""
        from workflow.lib.absence_detection import (
            classify_liftoff_genes,
            filter_missing_genes
        )

        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        missing = filter_missing_genes(results)

        # Should only include unmapped genes (missing status)
        assert len(missing) == 2
        assert all(g["status"] == "missing" for g in missing)

    def test_filter_missing_includes_uncertain_option(
        self, liftoff_gff_file: Path, unmapped_file: Path
    ):
        """Test filtering missing genes with option to include uncertain."""
        from workflow.lib.absence_detection import (
            classify_liftoff_genes,
            filter_missing_genes
        )

        results = classify_liftoff_genes(
            lifted_gff=liftoff_gff_file,
            unmapped_file=unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        # Include uncertain genes in missing list
        missing = filter_missing_genes(results, include_uncertain=True)

        # Should include both missing and uncertain
        assert len(missing) == 4  # 2 missing + 2 uncertain
        statuses = {g["status"] for g in missing}
        assert statuses == {"missing", "uncertain"}


# =============================================================================
# Test GFF with Different Attribute Formats
# =============================================================================

class TestDifferentGFFFormats:
    """Tests for handling different GFF attribute formats."""

    def test_liftoff_percentage_coverage(self, tmp_path: Path, empty_unmapped_file: Path):
        """Test handling coverage as percentage (0-100) vs decimal (0-1)."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        # Some tools output coverage as percentage
        gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene001;coverage=95
"""
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(gff_content)

        results = classify_liftoff_genes(
            lifted_gff=gff_file,
            unmapped_file=empty_unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        # Coverage > 1 should be normalized to 0-1 range
        gene001 = results[0]
        assert gene001["coverage"] == 0.95  # Normalized from 95

    def test_gff_with_extra_attributes(self, tmp_path: Path, empty_unmapped_file: Path):
        """Test GFF with additional attributes doesn't break parsing."""
        from workflow.lib.absence_detection import classify_liftoff_genes

        gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene001;Name=BRCA1;coverage=0.95;extra_attr=value;another=123
"""
        gff_file = tmp_path / "test.gff3"
        gff_file.write_text(gff_content)

        results = classify_liftoff_genes(
            lifted_gff=gff_file,
            unmapped_file=empty_unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        assert len(results) == 1
        assert results[0]["coverage"] == 0.95


# =============================================================================
# Test Gzip Compressed Input
# =============================================================================

class TestCompressedInput:
    """Tests for handling gzip-compressed input files."""

    def test_gzip_gff_file(self, liftoff_gff_content: str, tmp_path: Path, empty_unmapped_file: Path):
        """Test reading gzip-compressed GFF file."""
        import gzip
        from workflow.lib.absence_detection import classify_liftoff_genes

        gff_file = tmp_path / "lifted_annotation.gff3.gz"
        with gzip.open(gff_file, "wt", encoding="utf-8") as f:
            f.write(liftoff_gff_content)

        results = classify_liftoff_genes(
            lifted_gff=gff_file,
            unmapped_file=empty_unmapped_file,
            reference_species="human",
            target_species="mouse_lemur"
        )

        assert len(results) == 3


# =============================================================================
# Test Read Stats from TSV
# =============================================================================

class TestReadClassificationStatsFromTsv:
    """Tests for reading classification stats from TSV files."""

    def test_read_stats_from_tsv(self, tmp_path: Path):
        """Test reading stats from a classification TSV file."""
        from workflow.lib.absence_detection import read_classification_stats_from_tsv

        # Create a test TSV
        tsv_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene001\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.95\tNM_001
gene002\tTP53\thuman\tmouse_lemur\tuncertain\t0.42\t0.42\tNM_002
gene003\tFOXP2\thuman\tmouse_lemur\tmissing\t0.0\t0.0\tNM_003
gene004\tNOTCH1\thuman\tmouse_lemur\tmissing\t0.0\t0.0\tNM_004
"""
        tsv_file = tmp_path / "gene_classification.tsv"
        tsv_file.write_text(tsv_content)

        stats = read_classification_stats_from_tsv(tsv_file)

        assert stats["total_genes"] == 4
        assert stats["present"] == 1
        assert stats["uncertain"] == 1
        assert stats["missing"] == 2
        assert stats["present_rate"] == 0.25
        assert stats["missing_rate"] == 0.5

    def test_read_stats_from_missing_file(self, tmp_path: Path):
        """Test reading stats from non-existent file returns zeros."""
        from workflow.lib.absence_detection import read_classification_stats_from_tsv

        stats = read_classification_stats_from_tsv(tmp_path / "nonexistent.tsv")

        assert stats["total_genes"] == 0
        assert stats["present"] == 0
        assert stats["uncertain"] == 0
        assert stats["missing"] == 0
        assert stats["present_rate"] == 0.0
        assert stats["missing_rate"] == 0.0
