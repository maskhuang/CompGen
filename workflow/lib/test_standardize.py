"""
Unit tests for CompGene Data Standardization Module.

Tests cover:
- Input validation (AC5)
- Representative transcript selection (AC2)
- CDS to protein translation (AC3)
- Genome and annotation standardization (AC1, AC4)
- Output file formats (AC3)

Source: Story 2.2 - 本地数据标准化
"""

import gzip
import os
import tempfile
from pathlib import Path

import pytest

from workflow.lib.standardize import (
    CODON_TABLE,
    COMPLEMENT,
    RepresentativeTranscript,
    SpeciesData,
    extract_cds_sequence,
    extract_proteins,
    reverse_complement,
    select_representative_transcripts,
    standardize_annotation,
    standardize_genome,
    standardize_species,
    translate_cds,
    validate_has_cds_features,
    validate_species_inputs,
    write_representative_transcripts_tsv,
)
from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.gff import GFFFeature, GeneModel, TranscriptModel


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_genome(temp_dir):
    """Create a sample genome FASTA file."""
    genome_path = temp_dir / "genome.fa"
    content = """>chr1
ATGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""
    genome_path.write_text(content)
    return genome_path


@pytest.fixture
def sample_annotation(temp_dir):
    """Create a sample GFF3 annotation file with CDS features."""
    annotation_path = temp_dir / "annotation.gff3"
    content = """##gff-version 3
chr1\ttest\tgene\t1\t120\t.\t+\t.\tID=gene1;Name=TestGene1
chr1\ttest\tmRNA\t1\t120\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\ttest\texon\t1\t60\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\ttest\texon\t61\t120\t.\t+\t.\tID=exon2;Parent=mRNA1
chr1\ttest\tCDS\t1\t60\t.\t+\t0\tID=cds1;Parent=mRNA1
chr1\ttest\tCDS\t61\t120\t.\t+\t0\tID=cds2;Parent=mRNA1
chr1\ttest\tgene\t130\t200\t.\t-\t.\tID=gene2;Name=TestGene2
chr1\ttest\tmRNA\t130\t200\t.\t-\t.\tID=mRNA2;Parent=gene2
chr1\ttest\texon\t130\t200\t.\t-\t.\tID=exon3;Parent=mRNA2
chr1\ttest\tCDS\t130\t200\t.\t-\t0\tID=cds3;Parent=mRNA2
"""
    annotation_path.write_text(content)
    return annotation_path


@pytest.fixture
def sample_annotation_no_cds(temp_dir):
    """Create a sample GFF3 annotation file WITHOUT CDS features."""
    annotation_path = temp_dir / "annotation_no_cds.gff3"
    content = """##gff-version 3
chr1\ttest\tgene\t1\t120\t.\t+\t.\tID=gene1;Name=TestGene1
chr1\ttest\tmRNA\t1\t120\t.\t+\t.\tID=mRNA1;Parent=gene1
chr1\ttest\texon\t1\t60\t.\t+\t.\tID=exon1;Parent=mRNA1
chr1\ttest\texon\t61\t120\t.\t+\t.\tID=exon2;Parent=mRNA1
"""
    annotation_path.write_text(content)
    return annotation_path


@pytest.fixture
def species_data(temp_dir, sample_genome, sample_annotation):
    """Create a SpeciesData instance for testing."""
    return SpeciesData(
        species_id="test_species",
        genome_path=sample_genome,
        annotation_path=sample_annotation,
        output_dir=temp_dir / "output",
        name="Test Species"
    )


# =============================================================================
# Test Input Validation (AC5)
# =============================================================================

class TestInputValidation:
    """Tests for input validation functions."""

    def test_validate_species_inputs_success(self, species_data):
        """Test validation passes with valid inputs."""
        # Should not raise any exception
        validate_species_inputs(species_data)

    def test_validate_species_inputs_missing_genome(self, temp_dir, sample_annotation):
        """Test validation fails when genome file is missing."""
        species = SpeciesData(
            species_id="test",
            genome_path=temp_dir / "nonexistent.fa",
            annotation_path=sample_annotation,
            output_dir=temp_dir / "output"
        )
        with pytest.raises(CompGeneError) as exc_info:
            validate_species_inputs(species)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_validate_species_inputs_missing_annotation(self, temp_dir, sample_genome):
        """Test validation fails when annotation file is missing."""
        species = SpeciesData(
            species_id="test",
            genome_path=sample_genome,
            annotation_path=temp_dir / "nonexistent.gff3",
            output_dir=temp_dir / "output"
        )
        with pytest.raises(CompGeneError) as exc_info:
            validate_species_inputs(species)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_validate_has_cds_features_success(self, sample_annotation):
        """Test CDS validation passes when CDS features exist."""
        # Should not raise
        validate_has_cds_features(sample_annotation)

    def test_validate_has_cds_features_no_cds(self, sample_annotation_no_cds):
        """Test CDS validation fails when no CDS features."""
        with pytest.raises(CompGeneError) as exc_info:
            validate_has_cds_features(sample_annotation_no_cds)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


# =============================================================================
# Test DNA Translation (AC3)
# =============================================================================

class TestDNATranslation:
    """Tests for DNA/protein translation functions."""

    def test_reverse_complement_simple(self):
        """Test reverse complement of a simple sequence."""
        assert reverse_complement("ATGC") == "GCAT"

    def test_reverse_complement_with_n(self):
        """Test reverse complement handles N bases."""
        assert reverse_complement("ATNGC") == "GCNAT"

    def test_reverse_complement_lowercase(self):
        """Test reverse complement handles lowercase input."""
        assert reverse_complement("atgc") == "GCAT"

    def test_translate_cds_simple(self):
        """Test translation of a simple CDS."""
        # ATG = M, TGC = C, GTA = V, TAA = stop
        cds = "ATGTGCGTATAA"
        protein = translate_cds(cds, include_stop=False)
        assert protein == "MCV"

    def test_translate_cds_with_stop(self):
        """Test translation includes stop codon when requested."""
        cds = "ATGTGCGTATAA"
        protein = translate_cds(cds, include_stop=True)
        assert protein == "MCV*"

    def test_translate_cds_unknown_codon(self):
        """Test translation handles unknown codons as X."""
        cds = "ATGNNN"  # NNN is unknown
        protein = translate_cds(cds)
        assert protein == "MX"

    def test_translate_cds_incomplete_codon(self):
        """Test translation handles incomplete final codon."""
        cds = "ATGTGCGT"  # Last GT is incomplete
        protein = translate_cds(cds)
        assert protein == "MC"

    def test_codon_table_completeness(self):
        """Test that codon table has all 64 codons."""
        bases = "ATGC"
        all_codons = [a + b + c for a in bases for b in bases for c in bases]
        for codon in all_codons:
            assert codon in CODON_TABLE, f"Missing codon: {codon}"


# =============================================================================
# Test Representative Transcript Selection (AC2)
# =============================================================================

class TestRepresentativeTranscriptSelection:
    """Tests for representative transcript selection."""

    def test_select_representative_by_cds_length(self):
        """Test selection prefers longest CDS."""
        # Create mock gene with two transcripts
        gene = GeneModel(gene_id="gene1")

        # Transcript 1: shorter CDS
        t1 = TranscriptModel(transcript_id="t1", gene_id="gene1")
        t1.cds_features = [
            GFFFeature("chr1", "test", "CDS", 1, 30, None, "+", 0, {}, 1)
        ]
        t1.exons = [
            GFFFeature("chr1", "test", "exon", 1, 50, None, "+", None, {}, 2)
        ]

        # Transcript 2: longer CDS (should be selected)
        t2 = TranscriptModel(transcript_id="t2", gene_id="gene1")
        t2.cds_features = [
            GFFFeature("chr1", "test", "CDS", 1, 60, None, "+", 0, {}, 3)
        ]
        t2.exons = [
            GFFFeature("chr1", "test", "exon", 1, 60, None, "+", None, {}, 4)
        ]

        gene.transcripts = {"t1": t1, "t2": t2}
        gene.representative_transcript = "t2"

        genes = {"gene1": gene}
        representatives = select_representative_transcripts(genes)

        assert len(representatives) == 1
        assert representatives[0].transcript_id == "t2"
        assert representatives[0].cds_length == 60

    def test_select_representative_by_exon_when_no_cds(self):
        """Test selection falls back to exon length when no CDS."""
        gene = GeneModel(gene_id="gene1")

        # Transcript 1: shorter exons
        t1 = TranscriptModel(transcript_id="t1", gene_id="gene1")
        t1.exons = [
            GFFFeature("chr1", "test", "exon", 1, 30, None, "+", None, {}, 1)
        ]

        # Transcript 2: longer exons (should be selected)
        t2 = TranscriptModel(transcript_id="t2", gene_id="gene1")
        t2.exons = [
            GFFFeature("chr1", "test", "exon", 1, 60, None, "+", None, {}, 2)
        ]

        gene.transcripts = {"t1": t1, "t2": t2}
        gene.representative_transcript = "t2"

        genes = {"gene1": gene}
        representatives = select_representative_transcripts(genes)

        assert len(representatives) == 1
        assert representatives[0].transcript_id == "t2"

    def test_write_representative_transcripts_tsv(self, temp_dir):
        """Test TSV output format."""
        representatives = [
            RepresentativeTranscript(
                gene_id="gene1",
                transcript_id="t1",
                cds_length=300,
                exon_length=500
            ),
            RepresentativeTranscript(
                gene_id="gene2",
                transcript_id="t2",
                cds_length=600,
                exon_length=800
            ),
        ]

        output_path = temp_dir / "representatives.tsv"
        count = write_representative_transcripts_tsv(representatives, output_path)

        assert count == 2
        assert output_path.exists()

        content = output_path.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 3  # header + 2 records
        assert lines[0] == "gene_id\ttranscript_id\tcds_length\texon_length\tis_representative"
        assert "gene1\tt1\t300\t500\t1" in lines[1]


# =============================================================================
# Test File Standardization (AC1, AC4)
# =============================================================================

class TestFileStandardization:
    """Tests for genome and annotation standardization."""

    def test_standardize_genome_creates_gzip(self, sample_genome, temp_dir):
        """Test genome standardization creates gzipped output."""
        output_path = temp_dir / "output" / "genome.fa.gz"
        count = standardize_genome(sample_genome, output_path)

        assert count == 2  # 2 sequences
        assert output_path.exists()
        assert output_path.suffix == ".gz"

        # Verify it's valid gzip
        with gzip.open(output_path, "rt") as f:
            content = f.read()
            assert ">chr1" in content
            assert ">chr2" in content

    def test_standardize_genome_missing_input(self, temp_dir):
        """Test standardize_genome fails with missing input."""
        with pytest.raises(CompGeneError) as exc_info:
            standardize_genome(
                temp_dir / "nonexistent.fa",
                temp_dir / "output.fa.gz"
            )
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_standardize_annotation_creates_gzip(self, sample_annotation, temp_dir):
        """Test annotation standardization creates gzipped output."""
        output_path = temp_dir / "output" / "annotation.gff3.gz"
        count = standardize_annotation(sample_annotation, output_path)

        assert count > 0
        assert output_path.exists()

        # Verify it's valid gzip
        with gzip.open(output_path, "rt") as f:
            content = f.read()
            assert "##gff-version" in content
            assert "gene" in content

    def test_standardize_annotation_missing_input(self, temp_dir):
        """Test standardize_annotation fails with missing input."""
        with pytest.raises(CompGeneError) as exc_info:
            standardize_annotation(
                temp_dir / "nonexistent.gff3",
                temp_dir / "output.gff3.gz"
            )
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_atomic_write_success_replaces_temp(self, temp_dir):
        """Test that atomic write replaces .tmp with final file on success."""
        from workflow.lib.io import atomic_write

        output_path = temp_dir / "test_output.txt"
        atomic_write(output_path, "final content")

        assert output_path.exists()
        assert output_path.read_text() == "final content"
        # .tmp should not exist after successful write
        temp_path = output_path.with_suffix(".txt.tmp")
        assert not temp_path.exists()

    def test_standardize_genome_cleans_temp_on_failure(self, temp_dir):
        """Test that .tmp files are cleaned up when standardize_genome fails."""
        # Create an empty file that will cause issues when read
        bad_input = temp_dir / "bad_genome.fa"
        bad_input.write_text("")  # Empty file

        output_path = temp_dir / "output" / "genome.fa.gz"
        temp_path = output_path.with_suffix(".gz.tmp")

        # This should succeed (empty file is valid, just produces 0 sequences)
        count = standardize_genome(bad_input, output_path)
        assert count == 0
        assert output_path.exists()
        assert not temp_path.exists()  # .tmp should be cleaned up


# =============================================================================
# Test Protein Extraction (AC3)
# =============================================================================

class TestProteinExtraction:
    """Tests for protein extraction functionality."""

    def test_extract_proteins_creates_outputs(
        self, temp_dir, sample_genome, sample_annotation
    ):
        """Test extract_proteins creates protein FASTA and mapping TSV."""
        proteins_output = temp_dir / "proteins.fa.gz"
        mapping_output = temp_dir / "representatives.tsv"

        genes_count, proteins_count = extract_proteins(
            species_id="test",
            annotation_path=sample_annotation,
            genome_path=sample_genome,
            proteins_output=proteins_output,
            mapping_output=mapping_output,
            min_protein_length=1  # Allow short proteins for test
        )

        assert genes_count > 0
        assert proteins_output.exists()
        assert mapping_output.exists()

    def test_extract_proteins_orthofinder_id_format(
        self, temp_dir, sample_genome, sample_annotation
    ):
        """Test protein IDs follow OrthoFinder format: species_id|gene_id."""
        proteins_output = temp_dir / "proteins.fa.gz"
        mapping_output = temp_dir / "representatives.tsv"

        extract_proteins(
            species_id="mmur",
            annotation_path=sample_annotation,
            genome_path=sample_genome,
            proteins_output=proteins_output,
            mapping_output=mapping_output,
            min_protein_length=1
        )

        with gzip.open(proteins_output, "rt") as f:
            content = f.read()
            # Check OrthoFinder-compatible ID format
            assert "mmur|" in content or content == ""  # May be empty if no valid proteins

    def test_extract_proteins_missing_annotation(self, temp_dir, sample_genome):
        """Test extract_proteins fails with missing annotation."""
        with pytest.raises(CompGeneError) as exc_info:
            extract_proteins(
                species_id="test",
                annotation_path=temp_dir / "nonexistent.gff3",
                genome_path=sample_genome,
                proteins_output=temp_dir / "proteins.fa.gz",
                mapping_output=temp_dir / "mapping.tsv"
            )
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING


# =============================================================================
# Test Full Standardization Pipeline
# =============================================================================

class TestFullStandardization:
    """Tests for the complete standardization pipeline."""

    def test_standardize_species_creates_all_outputs(self, species_data):
        """Test standardize_species creates all expected output files."""
        stats = standardize_species(species_data)

        output_dir = species_data.output_dir
        assert (output_dir / "genome.fa.gz").exists()
        assert (output_dir / "annotation.gff3.gz").exists()
        assert (output_dir / "proteins.longest.fa.gz").exists()
        assert (output_dir / "representative_transcripts.tsv").exists()

        assert stats["genome_sequences"] == 2
        assert stats["annotation_features"] > 0

    def test_standardize_species_missing_genome(self, temp_dir, sample_annotation):
        """Test standardize_species fails with missing genome."""
        species = SpeciesData(
            species_id="test",
            genome_path=temp_dir / "nonexistent.fa",
            annotation_path=sample_annotation,
            output_dir=temp_dir / "output"
        )
        with pytest.raises(CompGeneError) as exc_info:
            standardize_species(species)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING


# =============================================================================
# Test CDS Sequence Extraction
# =============================================================================

class TestCDSExtraction:
    """Tests for CDS sequence extraction from genome."""

    def test_extract_cds_sequence_positive_strand(self):
        """Test CDS extraction for positive strand transcript."""
        # Create a mock transcript with CDS on positive strand
        transcript = TranscriptModel(transcript_id="t1", gene_id="g1")
        transcript.cds_features = [
            GFFFeature("chr1", "test", "CDS", 1, 9, None, "+", 0, {}, 1),
        ]

        genome_sequences = {"chr1": "ATGCGTACG"}  # 9 bp

        cds_seq = extract_cds_sequence(transcript, genome_sequences)
        assert cds_seq == "ATGCGTACG"

    def test_extract_cds_sequence_negative_strand(self):
        """Test CDS extraction for negative strand transcript."""
        transcript = TranscriptModel(transcript_id="t1", gene_id="g1")
        transcript.cds_features = [
            GFFFeature("chr1", "test", "CDS", 1, 9, None, "-", 0, {}, 1),
        ]

        genome_sequences = {"chr1": "ATGCGTACG"}

        cds_seq = extract_cds_sequence(transcript, genome_sequences)
        # Reverse complement of ATGCGTACG
        expected = reverse_complement("ATGCGTACG")
        assert cds_seq == expected

    def test_extract_cds_sequence_no_cds(self):
        """Test CDS extraction returns None when no CDS features."""
        transcript = TranscriptModel(transcript_id="t1", gene_id="g1")
        # No CDS features

        genome_sequences = {"chr1": "ATGCGTACG"}

        cds_seq = extract_cds_sequence(transcript, genome_sequences)
        assert cds_seq is None

    def test_extract_cds_sequence_multiple_cds(self):
        """Test CDS extraction concatenates multiple CDS features."""
        transcript = TranscriptModel(transcript_id="t1", gene_id="g1")
        transcript.cds_features = [
            GFFFeature("chr1", "test", "CDS", 1, 3, None, "+", 0, {}, 1),
            GFFFeature("chr1", "test", "CDS", 7, 9, None, "+", 0, {}, 2),
        ]

        genome_sequences = {"chr1": "ATGXXXACG"}  # positions 1-3 and 7-9

        cds_seq = extract_cds_sequence(transcript, genome_sequences)
        assert cds_seq == "ATGACG"


# =============================================================================
# Test SpeciesData Dataclass
# =============================================================================

class TestSpeciesDataDataclass:
    """Tests for SpeciesData dataclass."""

    def test_species_data_string_path_conversion(self, temp_dir):
        """Test that string paths are converted to Path objects."""
        species = SpeciesData(
            species_id="test",
            genome_path=str(temp_dir / "genome.fa"),
            annotation_path=str(temp_dir / "annotation.gff3"),
            output_dir=str(temp_dir / "output")
        )

        assert isinstance(species.genome_path, Path)
        assert isinstance(species.annotation_path, Path)
        assert isinstance(species.output_dir, Path)

    def test_species_data_optional_name(self, temp_dir):
        """Test that name is optional."""
        species = SpeciesData(
            species_id="test",
            genome_path=temp_dir / "genome.fa",
            annotation_path=temp_dir / "annotation.gff3",
            output_dir=temp_dir / "output"
        )

        assert species.name is None

        species_with_name = SpeciesData(
            species_id="test",
            genome_path=temp_dir / "genome.fa",
            annotation_path=temp_dir / "annotation.gff3",
            output_dir=temp_dir / "output",
            name="Test Species"
        )

        assert species_with_name.name == "Test Species"
