"""
Tests for CompGene Annotation Enhance Module.

Tests the GFF3 enhancement functionality that adds Liftoff classification
attributes to lifted annotation files.

Source: Story 5.3 - 迁移注释输出
"""

import csv
import pytest
from pathlib import Path
from tempfile import TemporaryDirectory

from workflow.lib.annotation_enhance import (
    load_classification_data,
    format_gff3_attributes,
    enhance_gff_line,
    filter_by_status,
    enhance_gff_with_liftoff_attrs,
)


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def sample_gff_content():
    """Sample GFF3 content for testing."""
    return """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95;sequence_ID=NM_001234
chr1\tliftoff\tmRNA\t1000\t2000\t.\t+\t.\tID=mrna1;Parent=gene1;coverage=0.95
chr1\tliftoff\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mrna1
chr1\tliftoff\tCDS\t1000\t1500\t.\t+\t0\tID=cds1;Parent=mrna1
chr1\tliftoff\tgene\t3000\t4000\t.\t-\t.\tID=gene2;coverage=0.45;sequence_ID=NM_005678
chr1\tliftoff\tmRNA\t3000\t4000\t.\t-\t.\tID=mrna2;Parent=gene2;coverage=0.45
chr1\tliftoff\tgene\t5000\t6000\t.\t+\t.\tID=gene3;coverage=0.0;sequence_ID=NM_009999
"""


@pytest.fixture
def sample_classification_content():
    """Sample gene_classification.tsv content."""
    return """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
gene2\tTP53\thuman\tmouse_lemur\tuncertain\t0.45\t0.42\tNM_005678
gene3\tFOXP2\thuman\tmouse_lemur\tmissing\t0.0\t0.0\tNM_009999
"""


@pytest.fixture
def temp_gff_file(sample_gff_content):
    """Create a temporary GFF3 file."""
    with TemporaryDirectory() as tmpdir:
        gff_path = Path(tmpdir) / "test.gff3"
        gff_path.write_text(sample_gff_content)
        yield gff_path


@pytest.fixture
def temp_classification_file(sample_classification_content):
    """Create a temporary classification TSV file."""
    with TemporaryDirectory() as tmpdir:
        tsv_path = Path(tmpdir) / "gene_classification.tsv"
        tsv_path.write_text(sample_classification_content)
        yield tsv_path


# =============================================================================
# Tests for load_classification_data
# =============================================================================

class TestLoadClassificationData:
    """Tests for load_classification_data function."""

    def test_load_valid_classification(self, temp_classification_file):
        """Test loading valid classification data."""
        data = load_classification_data(temp_classification_file)

        assert "gene1" in data
        assert data["gene1"]["status"] == "present"
        assert data["gene1"]["coverage"] == 0.95
        assert data["gene1"]["identity"] == 0.92
        assert data["gene1"]["source_gene"] == "NM_001234"

    def test_load_all_genes(self, temp_classification_file):
        """Test that all genes are loaded."""
        data = load_classification_data(temp_classification_file)

        assert len(data) == 3
        assert "gene1" in data
        assert "gene2" in data
        assert "gene3" in data

    def test_missing_file_raises_error(self):
        """Test that missing file raises CompGeneError."""
        from workflow.lib.errors import CompGeneError

        with pytest.raises(CompGeneError):
            load_classification_data(Path("/nonexistent/path.tsv"))


# =============================================================================
# Tests for format_gff3_attributes
# =============================================================================

class TestFormatGFF3Attributes:
    """Tests for format_gff3_attributes function."""

    def test_format_simple_attributes(self):
        """Test formatting simple attributes."""
        attrs = {"ID": "gene1", "Name": "BRCA1"}
        result = format_gff3_attributes(attrs)

        assert "ID=gene1" in result
        assert "Name=BRCA1" in result
        assert ";" in result

    def test_format_with_special_chars(self):
        """Test formatting attributes with special characters."""
        attrs = {"ID": "gene1", "Note": "Test;note"}
        result = format_gff3_attributes(attrs)

        # Special chars should be URL-encoded
        assert "ID=gene1" in result


# =============================================================================
# Tests for enhance_gff_line
# =============================================================================

class TestEnhanceGFFLine:
    """Tests for enhance_gff_line function."""

    def test_enhance_gene_line(self):
        """Test enhancing a gene line with liftoff attributes."""
        line = "chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95"
        classification = {
            "status": "present",
            "coverage": 0.95,
            "identity": 0.92,
            "source_gene": "NM_001234",
        }

        result = enhance_gff_line(line, classification, "human")

        assert "liftoff_coverage=0.95" in result
        assert "liftoff_identity=0.92" in result
        assert "liftoff_status=present" in result
        assert "source_gene=NM_001234" in result
        assert "reference_species=human" in result

    def test_enhance_preserves_original_attrs(self):
        """Test that original attributes are preserved."""
        line = "chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95;Name=BRCA1"
        classification = {
            "status": "present",
            "coverage": 0.95,
            "identity": 0.92,
            "source_gene": "NM_001234",
        }

        result = enhance_gff_line(line, classification, "human")

        assert "ID=gene1" in result
        assert "Name=BRCA1" in result
        assert "coverage=0.95" in result

    def test_enhance_child_feature(self):
        """Test enhancing child feature (mRNA/CDS) with parent's classification."""
        line = "chr1\tliftoff\tmRNA\t1000\t2000\t.\t+\t.\tID=mrna1;Parent=gene1"
        classification = {
            "status": "present",
            "coverage": 0.95,
            "identity": 0.92,
            "source_gene": "NM_001234",
        }

        result = enhance_gff_line(line, classification, "human")

        assert "liftoff_status=present" in result
        assert "Parent=gene1" in result


# =============================================================================
# Tests for filter_by_status
# =============================================================================

class TestFilterByStatus:
    """Tests for filter_by_status function."""

    def test_filter_present_only(self):
        """Test filtering to only 'present' genes."""
        statuses = ["present", "uncertain", "missing"]

        present = filter_by_status(statuses, include_uncertain=False)

        assert "present" in present
        assert "uncertain" not in present
        assert "missing" not in present

    def test_filter_include_uncertain(self):
        """Test filtering with include_uncertain=True."""
        statuses = ["present", "uncertain", "missing"]

        allowed = filter_by_status(statuses, include_uncertain=True)

        assert "present" in allowed
        assert "uncertain" in allowed
        assert "missing" not in allowed


# =============================================================================
# Tests for enhance_gff_with_liftoff_attrs
# =============================================================================

class TestEnhanceGFFWithLiftoffAttrs:
    """Tests for the main enhancement function."""

    def test_enhance_outputs_present_genes_only(self):
        """Test that by default only 'present' genes are output."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95;sequence_ID=NM_001234
chr1\tliftoff\tgene\t3000\t4000\t.\t-\t.\tID=gene2;coverage=0.45;sequence_ID=NM_005678
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
gene2\tTP53\thuman\tmouse_lemur\tuncertain\t0.45\t0.42\tNM_005678
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            stats = enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
                include_uncertain=False
            )

            # Read output and verify
            output_content = output_path.read_text()

            assert "gene1" in output_content
            assert "gene2" not in output_content
            assert stats["output_genes"] == 1
            assert stats["filtered_genes"] == 1

    def test_enhance_includes_uncertain_when_requested(self):
        """Test that uncertain genes are included when include_uncertain=True."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95;sequence_ID=NM_001234
chr1\tliftoff\tgene\t3000\t4000\t.\t-\t.\tID=gene2;coverage=0.45;sequence_ID=NM_005678
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
gene2\tTP53\thuman\tmouse_lemur\tuncertain\t0.45\t0.42\tNM_005678
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            stats = enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
                include_uncertain=True
            )

            output_content = output_path.read_text()

            assert "gene1" in output_content
            assert "gene2" in output_content
            assert stats["output_genes"] == 2

    def test_enhance_adds_liftoff_attributes(self):
        """Test that liftoff_* attributes are added."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            output_content = output_path.read_text()

            assert "liftoff_coverage=0.95" in output_content
            assert "liftoff_identity=0.92" in output_content
            assert "liftoff_status=present" in output_content
            assert "source_gene=NM_001234" in output_content
            assert "reference_species=human" in output_content

    def test_enhance_child_features_inherit_parent(self):
        """Test that child features inherit parent gene's classification."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
chr1\tliftoff\tmRNA\t1000\t2000\t.\t+\t.\tID=mrna1;Parent=gene1
chr1\tliftoff\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=mrna1
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            output_content = output_path.read_text()

            # All features should have liftoff attributes
            lines = [l for l in output_content.strip().split("\n") if not l.startswith("#")]
            for line in lines:
                assert "liftoff_status=present" in line

    def test_enhance_returns_statistics(self):
        """Test that statistics are returned correctly."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
chr1\tliftoff\tgene\t3000\t4000\t.\t-\t.\tID=gene2;coverage=0.45
chr1\tliftoff\tgene\t5000\t6000\t.\t+\t.\tID=gene3;coverage=0.0
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
gene2\tTP53\thuman\tmouse_lemur\tuncertain\t0.45\t0.42\tNM_005678
gene3\tFOXP2\thuman\tmouse_lemur\tmissing\t0.0\t0.0\tNM_009999
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            stats = enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
                include_uncertain=False
            )

            assert stats["input_genes"] == 3
            assert stats["output_genes"] == 1
            assert stats["filtered_genes"] == 2

    def test_enhance_empty_gff(self):
        """Test handling of empty GFF file."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            stats = enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            assert stats["input_genes"] == 0
            assert stats["output_genes"] == 0

    def test_enhance_gene_not_in_classification(self):
        """Test handling of genes not found in classification data."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
chr1\tliftoff\tgene\t3000\t4000\t.\t-\t.\tID=unknown_gene;coverage=0.80
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            stats = enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            # unknown_gene should be filtered (no classification = filtered)
            output_content = output_path.read_text()
            assert "gene1" in output_content
            assert "unknown_gene" not in output_content
            assert stats["unclassified_genes"] == 1

    def test_enhance_atomic_write(self):
        """Test that output is written atomically."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            # Output should exist and temp file should not
            assert output_path.exists()
            temp_path = output_path.with_suffix(".gff3.tmp")
            assert not temp_path.exists()


# =============================================================================
# Tests for Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_percentage_coverage_normalization(self):
        """Test that percentage coverage values are normalized to 0.0-1.0 range."""
        with TemporaryDirectory() as tmpdir:
            # Coverage as percentage (95.0 instead of 0.95)
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t95.0\t92.0\tNM_001234
"""
            tsv_path = Path(tmpdir) / "classification.tsv"
            tsv_path.write_text(classification_content)

            data = load_classification_data(tsv_path)

            # Values should be normalized to 0.0-1.0 range
            assert "gene1" in data
            assert data["gene1"]["coverage"] == 0.95
            assert data["gene1"]["identity"] == 0.92

    def test_gff_with_extra_attributes(self):
        """Test GFF with many extra attributes."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95;Name=BRCA1;Note=Test;Dbxref=GeneID:123
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            output_content = output_path.read_text()

            # All original attributes should be preserved
            assert "ID=gene1" in output_content
            assert "Name=BRCA1" in output_content
            assert "Note=Test" in output_content
            assert "Dbxref=GeneID:123" in output_content
            # New attributes should be added
            assert "liftoff_status=present" in output_content

    def test_gzip_input_support(self):
        """Test that gzip-compressed GFF3 input is handled correctly."""
        import gzip

        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
"""
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234
"""
            gff_path = Path(tmpdir) / "input.gff3.gz"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            # Write gzip-compressed GFF
            with gzip.open(gff_path, "wt") as f:
                f.write(gff_content)
            class_path.write_text(classification_content)

            stats = enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            output_content = output_path.read_text()

            assert "liftoff_status=present" in output_content
            assert stats["input_genes"] == 1
            assert stats["output_genes"] == 1

    def test_special_characters_in_source_gene(self):
        """Test that special characters in source_gene are URL-encoded in GFF3 output."""
        with TemporaryDirectory() as tmpdir:
            gff_content = """##gff-version 3
chr1\tliftoff\tgene\t1000\t2000\t.\t+\t.\tID=gene1;coverage=0.95
"""
            # source_gene with special characters (semicolon, equals, comma)
            classification_content = """gene_id\tgene_name\treference_species\ttarget_species\tstatus\tcoverage\tidentity\tsource_gene
gene1\tBRCA1\thuman\tmouse_lemur\tpresent\t0.95\t0.92\tNM_001234;alt=NM_999
"""
            gff_path = Path(tmpdir) / "input.gff3"
            class_path = Path(tmpdir) / "classification.tsv"
            output_path = Path(tmpdir) / "output.gff3"

            gff_path.write_text(gff_content)
            class_path.write_text(classification_content)

            enhance_gff_with_liftoff_attrs(
                gff_path=gff_path,
                classification_path=class_path,
                output_path=output_path,
                reference_species="human",
            )

            output_content = output_path.read_text()

            # Semicolon should be URL-encoded as %3B
            assert "source_gene=NM_001234%3Balt%3DNM_999" in output_content
            # The line should still be valid (9 fields)
            for line in output_content.split("\n"):
                if line and not line.startswith("#"):
                    assert len(line.split("\t")) == 9
