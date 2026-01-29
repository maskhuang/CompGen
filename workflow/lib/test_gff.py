"""
Unit tests for the GFF3/GTF parser module.

Tests cover:
- GFF3 parsing and attribute extraction
- GTF parsing and attribute extraction
- Gene hierarchy construction
- Representative transcript selection
- Format detection
- Error handling for malformed input
- gzip compression support

Source: Story 2.1 - GFF/FASTA 解析库
"""

import gzip
import tempfile
from pathlib import Path

import pytest

from workflow.lib.errors import ErrorCode, CompGeneError
from workflow.lib.gff import (
    GFFFeature,
    GeneModel,
    TranscriptModel,
    parse_gff3,
    parse_gtf,
    parse_annotation,
    parse_gff3_attributes,
    parse_gtf_attributes,
    parse_gff_line,
    validate_gff_line,
    build_gene_hierarchy,
    get_representative_transcript,
    detect_gff_format,
    detect_format,
    validate_gff,
    validate_file,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def gff3_content() -> str:
    """Sample GFF3 content for testing."""
    return """##gff-version 3
##sequence-region chr1 1 10000
chr1\tNCBI\tgene\t1000\t3000\t.\t+\t.\tID=gene001;Name=BRCA1
chr1\tNCBI\tmRNA\t1000\t3000\t.\t+\t.\tID=transcript001;Parent=gene001
chr1\tNCBI\texon\t1000\t1200\t.\t+\t.\tID=exon001;Parent=transcript001
chr1\tNCBI\texon\t1500\t1800\t.\t+\t.\tID=exon002;Parent=transcript001
chr1\tNCBI\tCDS\t1050\t1200\t.\t+\t0\tID=cds001;Parent=transcript001
chr1\tNCBI\tCDS\t1500\t1750\t.\t+\t0\tID=cds002;Parent=transcript001
"""


@pytest.fixture
def gtf_content() -> str:
    """Sample GTF content for testing."""
    return """chr1\tNCBI\tgene\t1000\t3000\t.\t+\t.\tgene_id "gene001"; gene_name "BRCA1";
chr1\tNCBI\ttranscript\t1000\t3000\t.\t+\t.\tgene_id "gene001"; transcript_id "transcript001";
chr1\tNCBI\texon\t1000\t1200\t.\t+\t.\tgene_id "gene001"; transcript_id "transcript001";
chr1\tNCBI\tCDS\t1050\t1200\t.\t+\t0\tgene_id "gene001"; transcript_id "transcript001";
"""


@pytest.fixture
def gff3_file(gff3_content: str, tmp_path: Path) -> Path:
    """Create a temporary GFF3 file."""
    path = tmp_path / "test.gff3"
    path.write_text(gff3_content)
    return path


@pytest.fixture
def gtf_file(gtf_content: str, tmp_path: Path) -> Path:
    """Create a temporary GTF file."""
    path = tmp_path / "test.gtf"
    path.write_text(gtf_content)
    return path


@pytest.fixture
def gff3_gzip_file(gff3_content: str, tmp_path: Path) -> Path:
    """Create a temporary gzip-compressed GFF3 file."""
    path = tmp_path / "test.gff3.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(gff3_content)
    return path


# =============================================================================
# Test GFF3 Attribute Parsing
# =============================================================================

class TestGFF3AttributeParsing:
    """Tests for GFF3 attribute parsing."""

    def test_parse_simple_attributes(self):
        """Test parsing simple key=value attributes."""
        attrs = parse_gff3_attributes("ID=gene001;Name=BRCA1")
        assert attrs["ID"] == "gene001"
        assert attrs["Name"] == "BRCA1"

    def test_parse_empty_attributes(self):
        """Test parsing empty attributes."""
        attrs = parse_gff3_attributes("")
        assert attrs == {}

        attrs = parse_gff3_attributes(".")
        assert attrs == {}

    def test_parse_url_encoded_attributes(self):
        """Test URL-decoding in attributes."""
        attrs = parse_gff3_attributes("Note=Gene%3Bwith%3Dsemicolon")
        assert attrs["Note"] == "Gene;with=semicolon"

    def test_parse_multiple_values(self):
        """Test parsing attributes with multiple values."""
        attrs = parse_gff3_attributes("ID=gene001;Parent=parent1,parent2")
        assert attrs["Parent"] == "parent1,parent2"


# =============================================================================
# Test GTF Attribute Parsing
# =============================================================================

class TestGTFAttributeParsing:
    """Tests for GTF attribute parsing."""

    def test_parse_simple_attributes(self):
        """Test parsing simple GTF attributes."""
        attrs = parse_gtf_attributes('gene_id "gene001"; transcript_id "transcript001";')
        assert attrs["gene_id"] == "gene001"
        assert attrs["transcript_id"] == "transcript001"

    def test_parse_empty_attributes(self):
        """Test parsing empty attributes."""
        attrs = parse_gtf_attributes("")
        assert attrs == {}

    def test_parse_attribute_with_spaces(self):
        """Test parsing attributes with extra spaces."""
        attrs = parse_gtf_attributes('gene_id "gene001" ; gene_name "BRCA1" ;')
        assert attrs["gene_id"] == "gene001"
        assert attrs["gene_name"] == "BRCA1"


# =============================================================================
# Test Line Validation
# =============================================================================

class TestLineValidation:
    """Tests for GFF/GTF line validation."""

    def test_valid_gff_line(self):
        """Test validation of a valid GFF line."""
        line = "chr1\tNCBI\tgene\t1000\t2000\t.\t+\t.\tID=gene001"
        is_valid, error = validate_gff_line(line, 1)
        assert is_valid
        assert error is None

    def test_invalid_field_count(self):
        """Test detection of wrong field count."""
        line = "chr1\tNCBI\tgene\t1000\t2000\t.\t+"
        is_valid, error = validate_gff_line(line, 1)
        assert not is_valid
        assert "9 tab-separated fields" in error

    def test_invalid_coordinates(self):
        """Test detection of non-integer coordinates."""
        line = "chr1\tNCBI\tgene\tXXX\t2000\t.\t+\t.\tID=gene001"
        is_valid, error = validate_gff_line(line, 1)
        assert not is_valid
        assert "integers" in error

    def test_start_greater_than_end(self):
        """Test detection of start > end."""
        line = "chr1\tNCBI\tgene\t3000\t2000\t.\t+\t.\tID=gene001"
        is_valid, error = validate_gff_line(line, 1)
        assert not is_valid
        assert "cannot be greater than end" in error

    def test_invalid_strand(self):
        """Test detection of invalid strand."""
        line = "chr1\tNCBI\tgene\t1000\t2000\t.\tX\t.\tID=gene001"
        is_valid, error = validate_gff_line(line, 1)
        assert not is_valid
        assert "Invalid strand" in error

    def test_comment_line(self):
        """Test that comment lines are valid."""
        is_valid, error = validate_gff_line("# comment", 1)
        assert is_valid
        assert error is None

    def test_empty_line(self):
        """Test that empty lines are valid."""
        is_valid, error = validate_gff_line("", 1)
        assert is_valid


# =============================================================================
# Test GFF3 Parsing
# =============================================================================

class TestGFF3Parsing:
    """Tests for GFF3 file parsing."""

    def test_parse_gff3_file(self, gff3_file: Path):
        """Test parsing a GFF3 file."""
        features = list(parse_gff3(gff3_file))

        # Should have 6 features (gene, mRNA, 2 exons, 2 CDS)
        assert len(features) == 6

        # Check gene feature
        gene = features[0]
        assert gene.type == "gene"
        assert gene.id == "gene001"
        assert gene.start == 1000
        assert gene.end == 3000
        assert gene.strand == "+"

    def test_parse_gff3_gzip(self, gff3_gzip_file: Path):
        """Test parsing a gzip-compressed GFF3 file."""
        features = list(parse_gff3(gff3_gzip_file))
        assert len(features) == 6

    def test_parse_missing_file(self, tmp_path: Path):
        """Test error handling for missing file."""
        with pytest.raises(CompGeneError) as exc_info:
            list(parse_gff3(tmp_path / "nonexistent.gff3"))
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_parse_invalid_format(self, tmp_path: Path):
        """Test error handling for invalid format."""
        bad_file = tmp_path / "bad.gff3"
        bad_file.write_text("chr1\tNCBI\tgene\tinvalid\t2000\t.\t+\t.\tID=gene001")

        with pytest.raises(CompGeneError) as exc_info:
            list(parse_gff3(bad_file))
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


# =============================================================================
# Test GTF Parsing
# =============================================================================

class TestGTFParsing:
    """Tests for GTF file parsing."""

    def test_parse_gtf_file(self, gtf_file: Path):
        """Test parsing a GTF file."""
        features = list(parse_gtf(gtf_file))

        # Should have 4 features
        assert len(features) == 4

        # Check gene feature
        gene = features[0]
        assert gene.type == "gene"
        assert gene.attributes["gene_id"] == "gene001"
        assert gene.attributes["gene_name"] == "BRCA1"


# =============================================================================
# Test Gene Hierarchy
# =============================================================================

class TestGeneHierarchy:
    """Tests for gene hierarchy construction."""

    def test_build_hierarchy(self, gff3_file: Path):
        """Test building gene hierarchy from GFF3."""
        features = parse_gff3(gff3_file)
        genes = build_gene_hierarchy(features)

        assert "gene001" in genes
        gene = genes["gene001"]

        # Check gene structure
        assert gene.gene_id == "gene001"
        assert len(gene.transcripts) == 1
        assert "transcript001" in gene.transcripts

        # Check transcript structure
        transcript = gene.transcripts["transcript001"]
        assert transcript.transcript_id == "transcript001"
        assert len(transcript.exons) == 2
        assert len(transcript.cds_features) == 2

    def test_representative_transcript_selection(self):
        """Test selection of representative transcript."""
        gene = GeneModel(gene_id="gene001")

        # Add two transcripts with different CDS lengths
        t1 = TranscriptModel(transcript_id="t1", gene_id="gene001")
        t1.cds_features = [
            GFFFeature("chr1", ".", "CDS", 1, 100, None, "+", 0, {}, 1),
        ]

        t2 = TranscriptModel(transcript_id="t2", gene_id="gene001")
        t2.cds_features = [
            GFFFeature("chr1", ".", "CDS", 1, 200, None, "+", 0, {}, 1),
        ]

        gene.transcripts["t1"] = t1
        gene.transcripts["t2"] = t2

        # t2 should be representative (longer CDS)
        rep_id = get_representative_transcript(gene)
        assert rep_id == "t2"

    def test_representative_transcript_tiebreaker(self):
        """Test tiebreaker for transcripts with equal length."""
        gene = GeneModel(gene_id="gene001")

        # Add two transcripts with same CDS length
        t1 = TranscriptModel(transcript_id="t_b", gene_id="gene001")
        t1.cds_features = [
            GFFFeature("chr1", ".", "CDS", 1, 100, None, "+", 0, {}, 1),
        ]

        t2 = TranscriptModel(transcript_id="t_a", gene_id="gene001")
        t2.cds_features = [
            GFFFeature("chr1", ".", "CDS", 1, 100, None, "+", 0, {}, 1),
        ]

        gene.transcripts["t_b"] = t1
        gene.transcripts["t_a"] = t2

        # t_a should be representative (lexicographically smaller)
        rep_id = get_representative_transcript(gene)
        assert rep_id == "t_a"


# =============================================================================
# Test Format Detection
# =============================================================================

class TestFormatDetection:
    """Tests for format detection."""

    def test_detect_gff3_by_extension(self, tmp_path: Path):
        """Test GFF3 detection by extension."""
        path = tmp_path / "test.gff3"
        path.write_text("chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1")
        assert detect_gff_format(path) == "gff3"

    def test_detect_gtf_by_extension(self, tmp_path: Path):
        """Test GTF detection by extension."""
        path = tmp_path / "test.gtf"
        path.write_text('chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id "g1";')
        assert detect_gff_format(path) == "gtf"

    def test_detect_gff3_by_header(self, tmp_path: Path):
        """Test GFF3 detection by header."""
        path = tmp_path / "test.txt"
        path.write_text("##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1")
        assert detect_gff_format(path) == "gff3"

    def test_detect_gtf_by_content(self, tmp_path: Path):
        """Test GTF detection by content."""
        path = tmp_path / "test.txt"
        path.write_text('chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";')
        assert detect_gff_format(path) == "gtf"

    def test_detect_format_fasta(self, tmp_path: Path):
        """Test FASTA format detection."""
        path = tmp_path / "test.fa"
        path.write_text(">seq1\nATCG")
        assert detect_format(path) == "fasta"

    def test_detect_format_unknown(self, tmp_path: Path):
        """Test error on unknown format."""
        path = tmp_path / "test.txt"
        path.write_text("some random content\nthat is not a known format")
        with pytest.raises(CompGeneError) as exc_info:
            detect_format(path)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


# =============================================================================
# Test Validation
# =============================================================================

class TestValidation:
    """Tests for file validation."""

    def test_validate_gff_valid(self, gff3_file: Path):
        """Test validation of valid GFF3 file."""
        is_valid, error = validate_gff(gff3_file)
        assert is_valid
        assert error is None

    def test_validate_gff_invalid(self, tmp_path: Path):
        """Test validation of invalid GFF3 file."""
        bad_file = tmp_path / "bad.gff3"
        bad_file.write_text("chr1\tNCBI\tgene\tinvalid\t2000\t.\t+\t.\tID=gene001")

        is_valid, error = validate_gff(bad_file)
        assert not is_valid
        assert error is not None

    def test_validate_file_format_mismatch(self, gff3_file: Path):
        """Test validation with format mismatch."""
        is_valid, error = validate_file(gff3_file, expected_format="gtf")
        assert not is_valid
        assert "Expected gtf" in error


# =============================================================================
# Test Auto-detection Parsing
# =============================================================================

class TestAutoDetectionParsing:
    """Tests for auto-detection parsing."""

    def test_parse_annotation_gff3(self, gff3_file: Path):
        """Test parse_annotation with GFF3."""
        features = list(parse_annotation(gff3_file))
        assert len(features) == 6

    def test_parse_annotation_gtf(self, gtf_file: Path):
        """Test parse_annotation with GTF."""
        features = list(parse_annotation(gtf_file))
        assert len(features) == 4


# =============================================================================
# Test GFFFeature Properties
# =============================================================================

class TestGFFFeatureProperties:
    """Tests for GFFFeature dataclass properties."""

    def test_feature_length(self):
        """Test feature length calculation."""
        feature = GFFFeature(
            seqid="chr1", source=".", type="gene",
            start=1000, end=2000, score=None,
            strand="+", phase=None, attributes={}, line_num=1
        )
        assert feature.length == 1001  # 2000 - 1000 + 1

    def test_feature_id_gff3(self):
        """Test ID property for GFF3."""
        feature = GFFFeature(
            seqid="chr1", source=".", type="gene",
            start=1, end=100, score=None,
            strand="+", phase=None,
            attributes={"ID": "gene001", "Name": "BRCA1"},
            line_num=1
        )
        assert feature.id == "gene001"
        assert feature.name == "BRCA1"

    def test_feature_id_gtf(self):
        """Test ID property for GTF (transcript_id)."""
        feature = GFFFeature(
            seqid="chr1", source=".", type="transcript",
            start=1, end=100, score=None,
            strand="+", phase=None,
            attributes={"gene_id": "g1", "transcript_id": "t1"},
            line_num=1
        )
        assert feature.id == "t1"
        assert feature.parent == "g1"


# =============================================================================
# Test Actual Fixture Files
# =============================================================================

class TestFixtureFiles:
    """Tests that verify the actual fixture files work correctly."""

    @pytest.fixture
    def fixtures_dir(self) -> Path:
        """Get the path to the fixtures directory."""
        return Path(__file__).parent.parent.parent / "tests" / "fixtures"

    def test_mini_annotation_gff3_fixture(self, fixtures_dir: Path):
        """Test parsing the actual mini_annotation.gff3 fixture file."""
        path = fixtures_dir / "mini_annotation.gff3"
        assert path.exists(), f"Fixture file not found: {path}"

        features = list(parse_gff3(path))

        # Should have features for genes, mRNAs, exons, CDS
        assert len(features) > 0

        # Check for expected feature types
        types = {f.type for f in features}
        assert "gene" in types
        assert "mRNA" in types
        assert "exon" in types
        assert "CDS" in types

    def test_mini_annotation_gff3_hierarchy(self, fixtures_dir: Path):
        """Test building hierarchy from mini_annotation.gff3 fixture."""
        path = fixtures_dir / "mini_annotation.gff3"
        features = parse_gff3(path)
        genes = build_gene_hierarchy(features)

        # Should have 2 genes: gene001 (BRCA1) and gene002 (TP53)
        assert len(genes) == 2
        assert "gene001" in genes
        assert "gene002" in genes

        # gene001 should have 2 transcripts
        gene001 = genes["gene001"]
        assert len(gene001.transcripts) == 2
        assert "transcript001" in gene001.transcripts
        assert "transcript002" in gene001.transcripts

        # Verify transcript001 has exons and CDS
        t1 = gene001.transcripts["transcript001"]
        assert len(t1.exons) == 3
        assert len(t1.cds_features) == 3

        # Verify representative transcript is selected (should be transcript001 - longer CDS)
        assert gene001.representative_transcript == "transcript001"

    def test_mini_annotation_gtf_fixture(self, fixtures_dir: Path):
        """Test parsing the actual mini_annotation.gtf fixture file."""
        path = fixtures_dir / "mini_annotation.gtf"
        assert path.exists(), f"Fixture file not found: {path}"

        features = list(parse_gtf(path))
        assert len(features) > 0

        # Check GTF attributes are parsed
        gene_features = [f for f in features if f.type == "gene"]
        assert len(gene_features) > 0
        assert "gene_id" in gene_features[0].attributes

    def test_mini_annotation_gtf_hierarchy(self, fixtures_dir: Path):
        """Test building hierarchy from mini_annotation.gtf fixture (Issue #1 fix)."""
        path = fixtures_dir / "mini_annotation.gtf"
        features = parse_gtf(path)
        genes = build_gene_hierarchy(features)

        # Should have at least 1 gene
        assert len(genes) >= 1
        assert "gene001" in genes

        gene001 = genes["gene001"]
        assert len(gene001.transcripts) >= 1

        # Verify that exons and CDS are correctly assigned to transcripts (Issue #1)
        for transcript in gene001.transcripts.values():
            # GTF fixture has exons and CDS - they should be assigned
            assert len(transcript.exons) > 0 or len(transcript.cds_features) > 0, \
                f"Transcript {transcript.transcript_id} has no exons or CDS assigned"

    def test_gff3_fixture_validation(self, fixtures_dir: Path):
        """Test that mini_annotation.gff3 passes validation."""
        path = fixtures_dir / "mini_annotation.gff3"
        is_valid, error = validate_gff(path)
        assert is_valid, f"Fixture validation failed: {error}"

    def test_gtf_fixture_validation(self, fixtures_dir: Path):
        """Test that mini_annotation.gtf passes validation."""
        path = fixtures_dir / "mini_annotation.gtf"
        is_valid, error = validate_gff(path)
        assert is_valid, f"Fixture validation failed: {error}"


# =============================================================================
# Test Multi-Parent Features (Issue #5)
# =============================================================================

class TestMultiParentFeatures:
    """Tests for GFF3 multi-parent feature handling."""

    def test_multi_parent_exon(self, tmp_path: Path):
        """Test that exons with multiple parents are assigned to all transcripts."""
        gff_content = """##gff-version 3
chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene001
chr1\t.\tmRNA\t1\t500\t.\t+\t.\tID=t1;Parent=gene001
chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=t2;Parent=gene001
chr1\t.\texon\t100\t200\t.\t+\t.\tID=exon1;Parent=t1,t2
chr1\t.\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=t1,t2
"""
        path = tmp_path / "multi_parent.gff3"
        path.write_text(gff_content)

        features = parse_gff3(path)
        genes = build_gene_hierarchy(features)

        gene = genes["gene001"]
        assert len(gene.transcripts) == 2

        # Both transcripts should have the shared exon
        t1 = gene.transcripts["t1"]
        t2 = gene.transcripts["t2"]

        assert len(t1.exons) == 1
        assert len(t2.exons) == 1
        assert len(t1.cds_features) == 1
        assert len(t2.cds_features) == 1
