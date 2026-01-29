"""
Unit tests for bio_utils module.

Tests cover:
- Biopython availability detection
- Genetic code table translation
- Multi-format ID parsing
- ID consistency checking
- Audit summary generation

Source: Story 2.2b - Biopython 工具封装
"""

import gzip
import json
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from workflow.lib.bio_utils import (
    # Availability
    check_biopython_available,
    require_biopython,
    # Translation
    CODON_TABLE_NAMES,
    STANDARD_CODON_TABLE,
    translate_with_table,
    # ID parsing
    IDFormat,
    detect_id_format,
    parse_fasta_id,
    extract_ids,
    extract_ids_with_mapping,
    # Consistency checking
    ComparisonResult,
    compare_id_sets,
    extract_ids_from_tsv,
    check_orthofinder_consistency,
    check_eggnog_consistency,
    write_consistency_report,
    # Audit summary
    compute_length_histogram,
    detect_sequence_type,
    generate_audit_summary,
    write_summary_json,
)
from workflow.lib.errors import CompGeneError, ErrorCode


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def simple_fasta(temp_dir):
    """Create a simple FASTA file."""
    fasta_path = temp_dir / "test.fa"
    content = """>gene_001 description text
ATGCGATCGATCG
>gene_002 another description
ATGCGATCGATCGATCGATCG
>gene_003
ATGC
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def simple_fasta_gz(temp_dir):
    """Create a gzipped FASTA file."""
    fasta_path = temp_dir / "test.fa.gz"
    content = b""">gene_001 description
ATGCGATCGATCG
>gene_002
ATGCGATCGATCGATCGATCG
"""
    with gzip.open(fasta_path, 'wb') as f:
        f.write(content)
    return fasta_path


@pytest.fixture
def uniprot_fasta(temp_dir):
    """Create a UniProt-style FASTA file."""
    fasta_path = temp_dir / "uniprot.fa"
    content = """>sp|P12345|GENE1_HUMAN Some protein
MVLSPADKTN
>sp|Q67890|GENE2_MOUSE Another protein
MKTAYIAKQR
>tr|A0A123|GENE3_RAT Unreviewed
MSLRGKPGPM
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def ncbi_fasta(temp_dir):
    """Create an NCBI-style FASTA file."""
    fasta_path = temp_dir / "ncbi.fa"
    content = """>gi|123456|ref|NP_001234.1| protein 1
MVLSPADKTN
>ref|NP_005678.2| protein 2
MKTAYIAKQR
>gb|AAA12345.1| protein 3
MSLRGKPGPM
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def orthofinder_fasta(temp_dir):
    """Create an OrthoFinder-style FASTA file."""
    fasta_path = temp_dir / "orthofinder.fa"
    content = """>mmur|MMUR_00001
MVLSPADKTN
>mmur|MMUR_00002
MKTAYIAKQR
>lcat|LCAT_00001
MSLRGKPGPM
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def protein_fasta(temp_dir):
    """Create a protein FASTA file."""
    fasta_path = temp_dir / "proteins.fa"
    content = """>prot_001
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
>prot_002
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKGF
>prot_003
MSLRGKPGPMPPLAIDMLLNFVASFVKA
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def orthogroups_tsv(temp_dir):
    """Create an OrthoFinder Orthogroups.tsv file."""
    tsv_path = temp_dir / "Orthogroups.tsv"
    content = """Orthogroup\tmmur\tlcat
OG0000001\tmmur|MMUR_00001, mmur|MMUR_00002\tlcat|LCAT_00001
OG0000002\tmmur|MMUR_00003\t
OG0000003\t\tlcat|LCAT_00002, lcat|LCAT_00003
"""
    tsv_path.write_text(content)
    return tsv_path


@pytest.fixture
def eggnog_annotations(temp_dir):
    """Create an eggNOG-mapper output file."""
    tsv_path = temp_dir / "annotations.emapper"
    content = """# eggNOG-mapper v2.1.0
#query\tseed_ortholog\tevalue\tscore\ttaxonomic_scope
gene_001\tKOG1234\t1e-50\t250\tEukaryota
gene_002\tKOG5678\t1e-40\t200\tEukaryota
gene_004\tKOG9999\t1e-30\t150\tEukaryota
"""
    tsv_path.write_text(content)
    return tsv_path


# =============================================================================
# Test: Biopython Availability
# =============================================================================

class TestBiopythonAvailability:
    """Tests for Biopython availability checking."""

    def test_check_biopython_available_returns_bool(self):
        """check_biopython_available should return a boolean."""
        result = check_biopython_available()
        assert isinstance(result, bool)

    def test_require_biopython_when_unavailable(self):
        """require_biopython should raise when Biopython is not available."""
        with mock.patch('workflow.lib.bio_utils._BIOPYTHON_AVAILABLE', False):
            with mock.patch('workflow.lib.bio_utils.check_biopython_available', return_value=False):
                with pytest.raises(CompGeneError) as exc_info:
                    require_biopython("test_function")
                assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND
                assert "biopython" in str(exc_info.value).lower()


# =============================================================================
# Test: Genetic Code Translation
# =============================================================================

class TestTranslation:
    """Tests for genetic code table translation."""

    def test_standard_codon_table_completeness(self):
        """Standard codon table should have all 64 codons."""
        assert len(STANDARD_CODON_TABLE) == 64

    def test_codon_table_names_validity(self):
        """CODON_TABLE_NAMES should contain valid NCBI table IDs."""
        assert 1 in CODON_TABLE_NAMES
        assert CODON_TABLE_NAMES[1] == "Standard"
        assert 2 in CODON_TABLE_NAMES
        assert "Mitochondrial" in CODON_TABLE_NAMES[2]

    def test_translate_standard_sequence(self):
        """translate_with_table should translate standard sequences."""
        # ATG -> M (Methionine), TGG -> W (Tryptophan), TAA -> * (Stop)
        result = translate_with_table("ATGTGGTAA", table_id=1)
        assert result == "MW*" or result == "MW"  # Depends on to_stop

    def test_translate_with_to_stop(self):
        """translate_with_table with to_stop should stop at stop codon."""
        result = translate_with_table("ATGTGGTAAGGG", table_id=1, to_stop=True)
        assert result == "MW"

    def test_translate_handles_n_bases(self):
        """translate_with_table should handle N bases as X."""
        result = translate_with_table("ATGNNN", table_id=1)
        assert "X" in result

    def test_translate_invalid_table_raises(self):
        """translate_with_table should raise for invalid table ID."""
        with pytest.raises(CompGeneError) as exc_info:
            translate_with_table("ATGTGG", table_id=999)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT

    def test_translate_converts_rna_to_dna(self):
        """translate_with_table should handle RNA sequences (U -> T)."""
        result = translate_with_table("AUGUGGUAA", table_id=1, to_stop=True)
        assert result == "MW"


# =============================================================================
# Test: ID Format Detection
# =============================================================================

class TestIDFormatDetection:
    """Tests for FASTA ID format detection."""

    def test_detect_uniprot_sp(self):
        """detect_id_format should recognize UniProt sp format."""
        assert detect_id_format("sp|P12345|GENE_HUMAN") == IDFormat.UNIPROT

    def test_detect_uniprot_tr(self):
        """detect_id_format should recognize UniProt tr format."""
        assert detect_id_format("tr|A0A123|GENE_MOUSE") == IDFormat.UNIPROT

    def test_detect_ncbi_gi(self):
        """detect_id_format should recognize NCBI GI format."""
        assert detect_id_format("gi|123456|ref|NP_001234.1|") == IDFormat.NCBI_GI

    def test_detect_ncbi_ref(self):
        """detect_id_format should recognize NCBI RefSeq format."""
        assert detect_id_format("ref|NP_001234.1|") == IDFormat.NCBI_REF

    def test_detect_genbank(self):
        """detect_id_format should recognize GenBank format."""
        assert detect_id_format("gb|AAA12345.1|") == IDFormat.GENBANK

    def test_detect_orthofinder(self):
        """detect_id_format should recognize OrthoFinder format."""
        assert detect_id_format("species|gene_id") == IDFormat.ORTHOFINDER

    def test_detect_simple(self):
        """detect_id_format should recognize simple format."""
        assert detect_id_format("gene_001 description") == IDFormat.SIMPLE
        assert detect_id_format("gene_001") == IDFormat.SIMPLE

    def test_detect_strips_prefix(self):
        """detect_id_format should handle > prefix."""
        assert detect_id_format(">sp|P12345|GENE") == IDFormat.UNIPROT


# =============================================================================
# Test: ID Parsing
# =============================================================================

class TestIDParsing:
    """Tests for FASTA ID parsing."""

    def test_parse_uniprot_id(self):
        """parse_fasta_id should extract UniProt accession."""
        assert parse_fasta_id("sp|P12345|GENE_HUMAN") == "P12345"
        assert parse_fasta_id("tr|A0A123|GENE_MOUSE") == "A0A123"

    def test_parse_ncbi_gi_id(self):
        """parse_fasta_id should extract RefSeq ID from GI format."""
        assert parse_fasta_id("gi|123456|ref|NP_001234.1|") == "NP_001234.1"

    def test_parse_ncbi_ref_id(self):
        """parse_fasta_id should extract RefSeq ID."""
        assert parse_fasta_id("ref|NP_001234.1|") == "NP_001234.1"

    def test_parse_genbank_id(self):
        """parse_fasta_id should extract GenBank accession."""
        assert parse_fasta_id("gb|AAA12345.1|") == "AAA12345.1"

    def test_parse_orthofinder_id(self):
        """parse_fasta_id should extract gene ID from OrthoFinder format."""
        assert parse_fasta_id("mmur|MMUR_00001") == "MMUR_00001"
        assert parse_fasta_id("species|gene_id") == "gene_id"

    def test_parse_simple_id(self):
        """parse_fasta_id should extract first token for simple format."""
        assert parse_fasta_id("gene_001 description text") == "gene_001"
        assert parse_fasta_id("gene_001") == "gene_001"

    def test_parse_empty_header(self):
        """parse_fasta_id should handle empty headers."""
        assert parse_fasta_id("") == ""
        assert parse_fasta_id(">") == ""

    def test_parse_with_format_hint(self):
        """parse_fasta_id should use format hint when provided."""
        # Force UniProt parsing on ambiguous input
        result = parse_fasta_id("a|b|c", format_hint=IDFormat.UNIPROT)
        assert result == "b"


class TestExtractIDs:
    """Tests for ID extraction from FASTA files."""

    def test_extract_ids_simple(self, simple_fasta):
        """extract_ids should extract all IDs from simple FASTA."""
        ids = extract_ids(simple_fasta)
        assert ids == {"gene_001", "gene_002", "gene_003"}

    def test_extract_ids_uniprot(self, uniprot_fasta):
        """extract_ids should extract UniProt accessions."""
        ids = extract_ids(uniprot_fasta)
        assert ids == {"P12345", "Q67890", "A0A123"}

    def test_extract_ids_ncbi(self, ncbi_fasta):
        """extract_ids should extract NCBI IDs."""
        ids = extract_ids(ncbi_fasta)
        assert "NP_001234.1" in ids
        assert "NP_005678.2" in ids
        assert "AAA12345.1" in ids

    def test_extract_ids_orthofinder(self, orthofinder_fasta):
        """extract_ids should extract gene IDs from OrthoFinder format."""
        ids = extract_ids(orthofinder_fasta)
        assert ids == {"MMUR_00001", "MMUR_00002", "LCAT_00001"}

    def test_extract_ids_gzip(self, simple_fasta_gz):
        """extract_ids should handle gzipped files."""
        ids = extract_ids(simple_fasta_gz)
        assert "gene_001" in ids
        assert "gene_002" in ids

    def test_extract_ids_with_mapping(self, orthofinder_fasta):
        """extract_ids_with_mapping should return IDs and mapping."""
        ids, mapping = extract_ids_with_mapping(orthofinder_fasta)
        assert "MMUR_00001" in ids
        assert "mmur|MMUR_00001" in mapping
        assert mapping["mmur|MMUR_00001"] == "MMUR_00001"


# =============================================================================
# Test: ID Consistency Checking
# =============================================================================

class TestComparisonResult:
    """Tests for ComparisonResult dataclass."""

    def test_match_rate_calculation(self):
        """ComparisonResult should calculate match rate correctly."""
        result = ComparisonResult(
            source_a="a", source_b="b",
            total_a=100, total_b=90, common=80,
            only_in_a=set(range(20)), only_in_b=set(range(10))
        )
        assert result.match_rate == 0.8

    def test_match_rate_empty(self):
        """ComparisonResult should handle empty source_a."""
        result = ComparisonResult(
            source_a="a", source_b="b",
            total_a=0, total_b=10, common=0
        )
        assert result.match_rate == 0.0

    def test_is_consistent(self):
        """ComparisonResult.is_consistent should check threshold."""
        result = ComparisonResult(
            source_a="a", source_b="b",
            total_a=100, total_b=100, common=96
        )
        assert result.is_consistent(threshold=0.95)
        assert not result.is_consistent(threshold=0.97)

    def test_to_dict(self):
        """ComparisonResult.to_dict should serialize correctly."""
        result = ComparisonResult(
            source_a="file_a", source_b="file_b",
            total_a=10, total_b=10, common=8,
            only_in_a={"a1", "a2"}, only_in_b={"b1", "b2"}
        )
        d = result.to_dict()
        assert d["source_a"] == "file_a"
        assert d["comparison"]["common"] == 8
        assert len(d["only_in_a"]) == 2
        assert d["only_in_a_truncated"] is False
        assert d["only_in_b_truncated"] is False

    def test_to_dict_truncation(self):
        """ComparisonResult.to_dict should indicate truncation."""
        # Create large sets that will be truncated
        large_set_a = {f"id_{i}" for i in range(150)}
        large_set_b = {f"other_{i}" for i in range(50)}
        result = ComparisonResult(
            source_a="a", source_b="b",
            total_a=200, total_b=100, common=50,
            only_in_a=large_set_a, only_in_b=large_set_b
        )
        d = result.to_dict(max_items=100)
        assert len(d["only_in_a"]) == 100
        assert d["only_in_a_truncated"] is True
        assert len(d["only_in_b"]) == 50
        assert d["only_in_b_truncated"] is False


class TestCompareIDSets:
    """Tests for compare_id_sets function."""

    def test_compare_identical_sets(self):
        """compare_id_sets should handle identical sets."""
        ids = {"a", "b", "c"}
        result = compare_id_sets(ids, ids)
        assert result.common == 3
        assert result.match_rate == 1.0
        assert len(result.only_in_a) == 0
        assert len(result.only_in_b) == 0

    def test_compare_disjoint_sets(self):
        """compare_id_sets should handle disjoint sets."""
        result = compare_id_sets({"a", "b"}, {"c", "d"})
        assert result.common == 0
        assert result.match_rate == 0.0

    def test_compare_overlapping_sets(self):
        """compare_id_sets should handle overlapping sets."""
        result = compare_id_sets({"a", "b", "c"}, {"b", "c", "d"})
        assert result.common == 2
        assert result.only_in_a == {"a"}
        assert result.only_in_b == {"d"}


class TestExtractIDsFromTSV:
    """Tests for extract_ids_from_tsv function."""

    def test_extract_from_tsv(self, temp_dir):
        """extract_ids_from_tsv should extract IDs from column."""
        tsv_path = temp_dir / "test.tsv"
        tsv_path.write_text("header\nid_001\nid_002\nid_003\n")
        ids = extract_ids_from_tsv(tsv_path, id_column=0)
        assert ids == {"id_001", "id_002", "id_003"}

    def test_extract_from_tsv_comma_separated(self, temp_dir):
        """extract_ids_from_tsv should handle comma-separated IDs."""
        tsv_path = temp_dir / "test.tsv"
        tsv_path.write_text("header\nid_001, id_002\nid_003\n")
        ids = extract_ids_from_tsv(tsv_path, id_column=0)
        assert ids == {"id_001", "id_002", "id_003"}


class TestConsistencyChecking:
    """Tests for consistency checking functions."""

    def test_check_orthofinder_consistency(self, orthofinder_fasta, orthogroups_tsv):
        """check_orthofinder_consistency should compare FASTA to Orthogroups."""
        result = check_orthofinder_consistency(orthofinder_fasta, orthogroups_tsv)
        # Input has 3 IDs, orthogroups has 5 unique gene IDs
        assert result.total_a == 3
        assert result.common >= 2  # At least MMUR_00001, MMUR_00002

    def test_check_eggnog_consistency(self, simple_fasta, eggnog_annotations):
        """check_eggnog_consistency should compare FASTA to eggNOG output."""
        result = check_eggnog_consistency(simple_fasta, eggnog_annotations)
        # simple_fasta has gene_001, gene_002, gene_003
        # eggnog has gene_001, gene_002, gene_004
        assert result.total_a == 3
        assert result.common == 2  # gene_001, gene_002

    def test_write_consistency_report(self, temp_dir):
        """write_consistency_report should write valid JSON."""
        result = ComparisonResult(
            source_a="a", source_b="b",
            total_a=10, total_b=10, common=9,
            only_in_a={"x"}, only_in_b={"y"}
        )
        output_path = temp_dir / "report.json"
        write_consistency_report(result, output_path, threshold=0.95)

        assert output_path.exists()
        with open(output_path) as f:
            report = json.load(f)
        assert report["status"] == "FAIL"  # 0.9 < 0.95
        assert "generated_at" in report


# =============================================================================
# Test: Audit Summary Generation
# =============================================================================

class TestLengthHistogram:
    """Tests for length histogram computation."""

    def test_compute_histogram_default_bins(self):
        """compute_length_histogram should use default bins."""
        lengths = [50, 150, 250, 750, 1500, 3000, 7000, 15000]
        result = compute_length_histogram(lengths)
        assert "bins" in result
        assert "counts" in result
        assert sum(result["counts"]) == len(lengths)

    def test_compute_histogram_custom_bins(self):
        """compute_length_histogram should accept custom bins."""
        lengths = [10, 20, 30, 40, 50]
        result = compute_length_histogram(lengths, bins=[0, 25, 50])
        assert result["bins"] == [0, 25, 50]
        # [0,25): 10,20; [25,50): 30,40; [50,inf): 50
        assert result["counts"] == [2, 2, 1]

    def test_compute_histogram_empty(self):
        """compute_length_histogram should handle empty list."""
        result = compute_length_histogram([])
        assert sum(result["counts"]) == 0


class TestSequenceTypeDetection:
    """Tests for sequence type detection."""

    def test_detect_nucleotide(self):
        """detect_sequence_type should identify nucleotide sequences."""
        sequences = ["ATGCGATCG", "ATGCATGCAT", "NNNATGCNNN"]
        assert detect_sequence_type(sequences) == "nucleotide"

    def test_detect_protein(self):
        """detect_sequence_type should identify protein sequences."""
        sequences = ["MVLSPADKTN", "MKTAYIAKQR", "MSLRGKPGPM"]
        assert detect_sequence_type(sequences) == "protein"

    def test_detect_empty(self):
        """detect_sequence_type should handle empty list."""
        assert detect_sequence_type([]) == "unknown"


class TestAuditSummary:
    """Tests for audit summary generation."""

    def test_generate_summary_nucleotide(self, simple_fasta):
        """generate_audit_summary should generate stats for nucleotide FASTA."""
        summary = generate_audit_summary(simple_fasta)

        assert summary["sequence_type"] == "nucleotide"
        assert summary["statistics"]["sequence_count"] == 3
        assert summary["statistics"]["total_length"] == 13 + 21 + 4  # 38
        assert summary["statistics"]["min_length"] == 4
        assert summary["statistics"]["max_length"] == 21
        assert "gc_content" in summary["statistics"]
        assert "length_histogram" in summary["statistics"]

    def test_generate_summary_protein(self, protein_fasta):
        """generate_audit_summary should generate stats for protein FASTA."""
        summary = generate_audit_summary(protein_fasta)

        assert summary["sequence_type"] == "protein"
        assert summary["statistics"]["sequence_count"] == 3
        assert summary["statistics"]["gc_content"] is None  # Not applicable

    def test_generate_summary_id_format(self, uniprot_fasta):
        """generate_audit_summary should detect ID format."""
        summary = generate_audit_summary(uniprot_fasta)

        assert summary["id_format"]["detected"] == IDFormat.UNIPROT
        assert len(summary["id_format"]["sample_ids"]) <= 3

    def test_write_summary_json(self, simple_fasta, temp_dir):
        """write_summary_json should write valid JSON file."""
        output_path = write_summary_json(simple_fasta)

        assert output_path.exists()
        assert output_path.name == "test.summary.json"

        with open(output_path) as f:
            summary = json.load(f)
        assert "statistics" in summary
        assert "generated_at" in summary

    def test_write_summary_json_custom_path(self, simple_fasta, temp_dir):
        """write_summary_json should accept custom output path."""
        custom_path = temp_dir / "custom_summary.json"
        result_path = write_summary_json(simple_fasta, custom_path)

        assert result_path == custom_path
        assert custom_path.exists()

    def test_write_summary_json_gz_extension(self, simple_fasta_gz):
        """write_summary_json should handle .fa.gz extension."""
        output_path = write_summary_json(simple_fasta_gz)
        assert output_path.name == "test.summary.json"


# =============================================================================
# Test: Integration
# =============================================================================

class TestIntegration:
    """Integration tests combining multiple functions."""

    def test_full_consistency_workflow(self, orthofinder_fasta, orthogroups_tsv, temp_dir):
        """Test complete consistency checking workflow."""
        # 1. Extract IDs from input
        input_ids = extract_ids(orthofinder_fasta)
        assert len(input_ids) == 3

        # 2. Check consistency with OrthoFinder output
        result = check_orthofinder_consistency(orthofinder_fasta, orthogroups_tsv)

        # 3. Write report
        report_path = temp_dir / "consistency_report.json"
        write_consistency_report(result, report_path)

        # 4. Verify report
        with open(report_path) as f:
            report = json.load(f)
        assert "status" in report
        assert "comparison" in report

    def test_full_audit_workflow(self, protein_fasta, temp_dir):
        """Test complete audit summary workflow."""
        # 1. Generate summary
        summary = generate_audit_summary(protein_fasta)
        assert summary["sequence_type"] == "protein"

        # 2. Write to file
        output_path = temp_dir / "proteins.summary.json"
        write_summary_json(protein_fasta, output_path)

        # 3. Read and verify
        with open(output_path) as f:
            loaded = json.load(f)
        assert loaded["statistics"]["sequence_count"] == summary["statistics"]["sequence_count"]
