"""
Unit tests for CompGene Orthogroup Utilities.

Tests cover:
- Parsing OrthoFinder Orthogroups.tsv to long format
- Orthogroup statistics calculation
- Species overlap matrix generation
- Edge cases (empty files, single species, missing genes)

Source: Story 3.2 - Orthogroups 表生成
"""

import csv
from pathlib import Path

import pytest

from workflow.lib.orthogroup_utils import (
    __version__,
    parse_orthogroups_tsv,
    calculate_orthogroup_stats,
    calculate_species_overlap,
    write_orthogroups_long_format,
    write_orthogroup_stats,
    write_species_overlap,
)


# =============================================================================
# Fixtures
# =============================================================================

SAMPLE_ORTHOGROUPS_TSV = """\
Orthogroup\tspecies1\tspecies2\tspecies3
OG0000000\tgene1, gene2\tgene3\tgene4, gene5, gene6
OG0000001\tgene7\tgene8\t
OG0000002\t\tgene9\tgene10
"""

SINGLE_COPY_ORTHOGROUPS_TSV = """\
Orthogroup\tspeciesA\tspeciesB
OG0000000\tgeneA1\tgeneB1
OG0000001\tgeneA2\tgeneB2
OG0000002\tgeneA3\t
"""

EMPTY_ORTHOGROUPS_TSV = """\
Orthogroup\tspecies1\tspecies2
"""


@pytest.fixture
def sample_orthogroups_file(tmp_path: Path) -> Path:
    path = tmp_path / "Orthogroups.tsv"
    path.write_text(SAMPLE_ORTHOGROUPS_TSV)
    return path


@pytest.fixture
def single_copy_file(tmp_path: Path) -> Path:
    path = tmp_path / "Orthogroups.tsv"
    path.write_text(SINGLE_COPY_ORTHOGROUPS_TSV)
    return path


@pytest.fixture
def empty_file(tmp_path: Path) -> Path:
    path = tmp_path / "Orthogroups.tsv"
    path.write_text(EMPTY_ORTHOGROUPS_TSV)
    return path


# =============================================================================
# Tests: parse_orthogroups_tsv
# =============================================================================

class TestParseOrthogroupsTsv:
    """Tests for parsing OrthoFinder Orthogroups.tsv to long format."""

    def test_basic_parsing(self, sample_orthogroups_file: Path):
        """Correctly parses multi-species orthogroups with comma-separated genes."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)

        # OG0000000: 2 + 1 + 3 = 6 genes
        # OG0000001: 1 + 1 + 0 = 2 genes
        # OG0000002: 0 + 1 + 1 = 2 genes
        assert len(records) == 10

    def test_record_fields(self, sample_orthogroups_file: Path):
        """Each record has orthogroup_id, gene_id, species_id."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        for rec in records:
            assert "orthogroup_id" in rec
            assert "gene_id" in rec
            assert "species_id" in rec

    def test_gene_species_mapping(self, sample_orthogroups_file: Path):
        """Genes are correctly mapped to their species."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        by_gene = {r["gene_id"]: r for r in records}

        assert by_gene["gene1"]["species_id"] == "species1"
        assert by_gene["gene1"]["orthogroup_id"] == "OG0000000"
        assert by_gene["gene3"]["species_id"] == "species2"
        assert by_gene["gene4"]["species_id"] == "species3"

    def test_empty_cells_skipped(self, sample_orthogroups_file: Path):
        """Empty cells (species with no gene in OG) produce no records."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        og0001_species = {r["species_id"] for r in records if r["orthogroup_id"] == "OG0000001"}
        assert og0001_species == {"species1", "species2"}
        # species3 has empty cell in OG0000001

    def test_empty_file(self, empty_file: Path):
        """Empty file (header only) returns empty list."""
        records = parse_orthogroups_tsv(empty_file)
        assert records == []

    def test_single_copy(self, single_copy_file: Path):
        """Single-copy orthogroups parse correctly."""
        records = parse_orthogroups_tsv(single_copy_file)
        og0000 = [r for r in records if r["orthogroup_id"] == "OG0000000"]
        assert len(og0000) == 2
        assert {r["gene_id"] for r in og0000} == {"geneA1", "geneB1"}

    def test_whitespace_handling(self, tmp_path: Path):
        """Gene IDs with extra whitespace are stripped."""
        content = "Orthogroup\tsp1\tsp2\nOG0\t gene_a , gene_b \tgene_c\n"
        path = tmp_path / "ws.tsv"
        path.write_text(content)
        records = parse_orthogroups_tsv(path)
        gene_ids = {r["gene_id"] for r in records}
        assert gene_ids == {"gene_a", "gene_b", "gene_c"}

    def test_file_not_found(self, tmp_path: Path):
        """Missing file raises CompGeneError."""
        from workflow.lib.errors import CompGeneError
        with pytest.raises(CompGeneError):
            parse_orthogroups_tsv(tmp_path / "nonexistent.tsv")

    def test_species_names_from_header(self, sample_orthogroups_file: Path):
        """Species names extracted from header row."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        species = {r["species_id"] for r in records}
        assert species == {"species1", "species2", "species3"}


# =============================================================================
# Tests: calculate_orthogroup_stats
# =============================================================================

class TestCalculateOrthogroupStats:
    """Tests for orthogroup statistics calculation."""

    def test_basic_stats(self, sample_orthogroups_file: Path):
        """Correctly calculates gene_count, species_count, is_single_copy."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        stats = calculate_orthogroup_stats(records)

        og0 = {s["orthogroup_id"]: s for s in stats}

        # OG0000000: 6 genes, 3 species, not single copy
        assert og0["OG0000000"]["gene_count"] == 6
        assert og0["OG0000000"]["species_count"] == 3
        assert og0["OG0000000"]["is_single_copy"] is False

        # OG0000001: 2 genes, 2 species, single copy (1 per species)
        assert og0["OG0000001"]["gene_count"] == 2
        assert og0["OG0000001"]["species_count"] == 2
        assert og0["OG0000001"]["is_single_copy"] is True

        # OG0000002: 2 genes, 2 species, single copy
        assert og0["OG0000002"]["gene_count"] == 2
        assert og0["OG0000002"]["species_count"] == 2
        assert og0["OG0000002"]["is_single_copy"] is True

    def test_species_list(self, sample_orthogroups_file: Path):
        """Species list for each orthogroup is correct."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        stats = calculate_orthogroup_stats(records)
        og0 = {s["orthogroup_id"]: s for s in stats}

        assert set(og0["OG0000000"]["species_list"].split(",")) == {"species1", "species2", "species3"}
        assert set(og0["OG0000002"]["species_list"].split(",")) == {"species2", "species3"}

    def test_empty_records(self):
        """Empty records list produces empty stats."""
        stats = calculate_orthogroup_stats([])
        assert stats == []

    def test_single_copy_detection(self, single_copy_file: Path):
        """Correctly identifies single-copy orthogroups."""
        records = parse_orthogroups_tsv(single_copy_file)
        stats = calculate_orthogroup_stats(records)
        og0 = {s["orthogroup_id"]: s for s in stats}

        assert og0["OG0000000"]["is_single_copy"] is True
        assert og0["OG0000001"]["is_single_copy"] is True
        # OG0000002 has only 1 species → not single copy (needs >=2 species for true single-copy)
        assert og0["OG0000002"]["is_single_copy"] is False


# =============================================================================
# Tests: calculate_species_overlap
# =============================================================================

class TestCalculateSpeciesOverlap:
    """Tests for species overlap matrix calculation."""

    def test_basic_overlap(self, sample_orthogroups_file: Path):
        """Correctly computes shared orthogroup counts."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        matrix = calculate_species_overlap(records)

        # species1 participates in OG0000000, OG0000001 → 2 OGs
        assert matrix["species1"]["species1"] == 2
        # species2 participates in OG0000000, OG0000001, OG0000002 → 3 OGs
        assert matrix["species2"]["species2"] == 3
        # species3 participates in OG0000000, OG0000002 → 2 OGs
        assert matrix["species3"]["species3"] == 2

    def test_shared_counts(self, sample_orthogroups_file: Path):
        """Pairwise shared orthogroup counts are correct."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        matrix = calculate_species_overlap(records)

        # species1 & species2 share: OG0000000, OG0000001 → 2
        assert matrix["species1"]["species2"] == 2
        # species1 & species3 share: OG0000000 → 1
        assert matrix["species1"]["species3"] == 1
        # species2 & species3 share: OG0000000, OG0000002 → 2
        assert matrix["species2"]["species3"] == 2

    def test_symmetry(self, sample_orthogroups_file: Path):
        """Matrix is symmetric."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        matrix = calculate_species_overlap(records)
        species = sorted(matrix.keys())
        for s1 in species:
            for s2 in species:
                assert matrix[s1][s2] == matrix[s2][s1]

    def test_empty_records(self):
        """Empty records produce empty matrix."""
        matrix = calculate_species_overlap([])
        assert matrix == {}


# =============================================================================
# Tests: write functions
# =============================================================================

class TestWriteFunctions:
    """Tests for file output functions."""

    def test_write_orthogroups_long_format(self, sample_orthogroups_file: Path, tmp_path: Path):
        """Write long format TSV with correct columns."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        output = tmp_path / "orthogroups.tsv"
        write_orthogroups_long_format(records, output)

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header == ["orthogroup_id", "gene_id", "species_id"]
        assert len(lines) == 11  # 1 header + 10 data rows

    def test_write_orthogroup_stats(self, sample_orthogroups_file: Path, tmp_path: Path):
        """Write stats TSV with correct columns."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        stats = calculate_orthogroup_stats(records)
        output = tmp_path / "orthogroup_stats.tsv"
        write_orthogroup_stats(stats, output)

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header == ["orthogroup_id", "gene_count", "species_count", "is_single_copy", "species_list"]
        assert len(lines) == 4  # 1 header + 3 OGs

    def test_write_species_overlap(self, sample_orthogroups_file: Path, tmp_path: Path):
        """Write overlap matrix TSV."""
        records = parse_orthogroups_tsv(sample_orthogroups_file)
        matrix = calculate_species_overlap(records)
        output = tmp_path / "species_overlap.tsv"
        write_species_overlap(matrix, output)

        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 4  # 1 header + 3 species

    def test_write_species_overlap_empty(self, tmp_path: Path):
        """Empty matrix writes header-only file."""
        output = tmp_path / "empty_overlap.tsv"
        write_species_overlap({}, output)
        assert output.exists()
        assert output.read_text() == "species\n"

    def test_atomic_write(self, tmp_path: Path):
        """No .tmp file left behind after successful write."""
        records = [{"orthogroup_id": "OG0", "gene_id": "g1", "species_id": "sp1"}]
        output = tmp_path / "test.tsv"
        write_orthogroups_long_format(records, output)
        tmp_file = output.with_suffix(".tsv.tmp")
        assert not tmp_file.exists()
        assert output.exists()
