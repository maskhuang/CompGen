"""
Unit tests for presence_absence module.

Tests cover:
- Binary presence/absence matrix generation
- Copy number matrix generation
- Matrix TSV writing
- Species list extraction
- Edge cases (empty input, single species, single orthogroup)
- Matrix consistency (PA vs CN)

Source: Story 6A.1
"""

import pytest
from pathlib import Path

from workflow.lib.presence_absence import (
    build_presence_absence_matrix,
    build_copy_number_matrix,
    get_sorted_species,
    read_orthogroups_long_format,
    write_matrix_tsv,
)
from workflow.lib.errors import CompGeneError


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def multi_species_records():
    """Records with multiple species and multi-copy genes."""
    return [
        {"orthogroup_id": "OG0000000", "gene_id": "g1", "species_id": "species1"},
        {"orthogroup_id": "OG0000000", "gene_id": "g2", "species_id": "species1"},
        {"orthogroup_id": "OG0000000", "gene_id": "g3", "species_id": "species2"},
        {"orthogroup_id": "OG0000000", "gene_id": "g4", "species_id": "species3"},
        {"orthogroup_id": "OG0000001", "gene_id": "g5", "species_id": "species1"},
        {"orthogroup_id": "OG0000001", "gene_id": "g6", "species_id": "species2"},
        {"orthogroup_id": "OG0000002", "gene_id": "g7", "species_id": "species2"},
        {"orthogroup_id": "OG0000002", "gene_id": "g8", "species_id": "species3"},
    ]


@pytest.fixture
def single_species_records():
    """Records with only one species."""
    return [
        {"orthogroup_id": "OG0000000", "gene_id": "g1", "species_id": "speciesA"},
        {"orthogroup_id": "OG0000001", "gene_id": "g2", "species_id": "speciesA"},
        {"orthogroup_id": "OG0000001", "gene_id": "g3", "species_id": "speciesA"},
    ]


@pytest.fixture
def single_orthogroup_records():
    """Records with only one orthogroup."""
    return [
        {"orthogroup_id": "OG0000000", "gene_id": "g1", "species_id": "sp1"},
        {"orthogroup_id": "OG0000000", "gene_id": "g2", "species_id": "sp2"},
        {"orthogroup_id": "OG0000000", "gene_id": "g3", "species_id": "sp2"},
        {"orthogroup_id": "OG0000000", "gene_id": "g4", "species_id": "sp3"},
    ]


# =============================================================================
# Tests: get_sorted_species
# =============================================================================

class TestGetSortedSpecies:
    def test_multiple_species(self, multi_species_records):
        result = get_sorted_species(multi_species_records)
        assert result == ["species1", "species2", "species3"]

    def test_empty_records(self):
        result = get_sorted_species([])
        assert result == []

    def test_single_species(self, single_species_records):
        result = get_sorted_species(single_species_records)
        assert result == ["speciesA"]

    def test_alphabetical_order(self):
        records = [
            {"orthogroup_id": "OG0", "gene_id": "g1", "species_id": "zebra"},
            {"orthogroup_id": "OG0", "gene_id": "g2", "species_id": "alpha"},
            {"orthogroup_id": "OG0", "gene_id": "g3", "species_id": "mouse"},
        ]
        result = get_sorted_species(records)
        assert result == ["alpha", "mouse", "zebra"]


# =============================================================================
# Tests: build_presence_absence_matrix
# =============================================================================

class TestBuildPresenceAbsenceMatrix:
    def test_basic_binary(self, multi_species_records):
        matrix = build_presence_absence_matrix(multi_species_records)
        # OG0000000: species1=1, species2=1, species3=1
        assert matrix["OG0000000"]["species1"] == 1
        assert matrix["OG0000000"]["species2"] == 1
        assert matrix["OG0000000"]["species3"] == 1
        # OG0000001: species1=1, species2=1, species3=0
        assert matrix["OG0000001"]["species1"] == 1
        assert matrix["OG0000001"]["species2"] == 1
        assert matrix["OG0000001"]["species3"] == 0
        # OG0000002: species1=0, species2=1, species3=1
        assert matrix["OG0000002"]["species1"] == 0
        assert matrix["OG0000002"]["species2"] == 1
        assert matrix["OG0000002"]["species3"] == 1

    def test_multi_copy_still_binary(self, multi_species_records):
        """Multi-copy genes should still show 1, not count."""
        matrix = build_presence_absence_matrix(multi_species_records)
        # species1 has 2 genes in OG0000000 but PA should be 1
        assert matrix["OG0000000"]["species1"] == 1

    def test_empty_records(self):
        matrix = build_presence_absence_matrix([])
        assert matrix == {}

    def test_single_species(self, single_species_records):
        matrix = build_presence_absence_matrix(single_species_records)
        assert matrix["OG0000000"]["speciesA"] == 1
        assert matrix["OG0000001"]["speciesA"] == 1

    def test_single_orthogroup(self, single_orthogroup_records):
        matrix = build_presence_absence_matrix(single_orthogroup_records)
        assert len(matrix) == 1
        assert matrix["OG0000000"]["sp1"] == 1
        assert matrix["OG0000000"]["sp2"] == 1
        assert matrix["OG0000000"]["sp3"] == 1

    def test_all_species_present_in_every_og(self, multi_species_records):
        """All species columns present in every orthogroup row (0 for absent)."""
        matrix = build_presence_absence_matrix(multi_species_records)
        species = {"species1", "species2", "species3"}
        for og_id, row in matrix.items():
            assert set(row.keys()) == species, f"{og_id} missing species columns"


# =============================================================================
# Tests: build_copy_number_matrix
# =============================================================================

class TestBuildCopyNumberMatrix:
    def test_basic_counts(self, multi_species_records):
        matrix = build_copy_number_matrix(multi_species_records)
        # OG0000000: species1=2 (g1, g2), species2=1 (g3), species3=1 (g4)
        assert matrix["OG0000000"]["species1"] == 2
        assert matrix["OG0000000"]["species2"] == 1
        assert matrix["OG0000000"]["species3"] == 1
        # OG0000001: species1=1, species2=1, species3=0
        assert matrix["OG0000001"]["species1"] == 1
        assert matrix["OG0000001"]["species2"] == 1
        assert matrix["OG0000001"]["species3"] == 0
        # OG0000002: species1=0, species2=1, species3=1
        assert matrix["OG0000002"]["species1"] == 0
        assert matrix["OG0000002"]["species2"] == 1
        assert matrix["OG0000002"]["species3"] == 1

    def test_empty_records(self):
        matrix = build_copy_number_matrix([])
        assert matrix == {}

    def test_single_species_multi_copy(self, single_species_records):
        matrix = build_copy_number_matrix(single_species_records)
        assert matrix["OG0000000"]["speciesA"] == 1
        assert matrix["OG0000001"]["speciesA"] == 2  # g2, g3

    def test_all_species_present_in_every_og(self, multi_species_records):
        """All species columns present in every orthogroup row (0 for absent)."""
        matrix = build_copy_number_matrix(multi_species_records)
        species = {"species1", "species2", "species3"}
        for og_id, row in matrix.items():
            assert set(row.keys()) == species


# =============================================================================
# Tests: Matrix Consistency (PA vs CN)
# =============================================================================

class TestMatrixConsistency:
    def test_pa_cn_consistency(self, multi_species_records):
        """PA matrix 1 ↔ CN matrix >=1, PA matrix 0 ↔ CN matrix 0."""
        pa = build_presence_absence_matrix(multi_species_records)
        cn = build_copy_number_matrix(multi_species_records)
        assert set(pa.keys()) == set(cn.keys()), "Same orthogroups"
        for og_id in pa:
            assert set(pa[og_id].keys()) == set(cn[og_id].keys()), "Same species"
            for sp in pa[og_id]:
                if pa[og_id][sp] == 1:
                    assert cn[og_id][sp] >= 1
                else:
                    assert cn[og_id][sp] == 0

    def test_pa_cn_same_dimensions(self, multi_species_records):
        pa = build_presence_absence_matrix(multi_species_records)
        cn = build_copy_number_matrix(multi_species_records)
        assert len(pa) == len(cn)
        for og_id in pa:
            assert len(pa[og_id]) == len(cn[og_id])


# =============================================================================
# Tests: write_matrix_tsv
# =============================================================================

class TestWriteMatrixTsv:
    def test_basic_write(self, multi_species_records, tmp_path: Path):
        matrix = build_presence_absence_matrix(multi_species_records)
        species = get_sorted_species(multi_species_records)
        output = tmp_path / "pa.tsv"
        write_matrix_tsv(matrix, species, output)
        assert output.exists()
        lines = output.read_text().strip().split("\n")
        # Header + 3 OGs
        assert len(lines) == 4
        header = lines[0].split("\t")
        assert header == ["orthogroup_id", "species1", "species2", "species3"]
        # OG0000000 row
        row0 = lines[1].split("\t")
        assert row0[0] == "OG0000000"
        assert row0 == ["OG0000000", "1", "1", "1"]

    def test_copy_number_write(self, multi_species_records, tmp_path: Path):
        matrix = build_copy_number_matrix(multi_species_records)
        species = get_sorted_species(multi_species_records)
        output = tmp_path / "cn.tsv"
        write_matrix_tsv(matrix, species, output)
        lines = output.read_text().strip().split("\n")
        # OG0000000: species1=2, species2=1, species3=1
        row0 = lines[1].split("\t")
        assert row0 == ["OG0000000", "2", "1", "1"]

    def test_empty_matrix(self, tmp_path: Path):
        output = tmp_path / "empty.tsv"
        write_matrix_tsv({}, [], output)
        assert output.exists()
        assert output.read_text() == "orthogroup_id\n"

    def test_atomic_write(self, multi_species_records, tmp_path: Path):
        """No .tmp file should remain after successful write."""
        matrix = build_presence_absence_matrix(multi_species_records)
        species = get_sorted_species(multi_species_records)
        output = tmp_path / "pa.tsv"
        write_matrix_tsv(matrix, species, output)
        tmp_file = output.with_suffix(".tsv.tmp")
        assert not tmp_file.exists()
        assert output.exists()

    def test_creates_parent_dirs(self, multi_species_records, tmp_path: Path):
        matrix = build_presence_absence_matrix(multi_species_records)
        species = get_sorted_species(multi_species_records)
        output = tmp_path / "sub" / "dir" / "pa.tsv"
        write_matrix_tsv(matrix, species, output)
        assert output.exists()

    def test_sorted_orthogroup_rows(self, multi_species_records, tmp_path: Path):
        """Orthogroup rows should be sorted by orthogroup_id."""
        matrix = build_presence_absence_matrix(multi_species_records)
        species = get_sorted_species(multi_species_records)
        output = tmp_path / "pa.tsv"
        write_matrix_tsv(matrix, species, output)
        lines = output.read_text().strip().split("\n")
        og_ids = [line.split("\t")[0] for line in lines[1:]]
        assert og_ids == sorted(og_ids)


# =============================================================================
# Tests: read_orthogroups_long_format
# =============================================================================

class TestReadOrthogroupsLongFormat:
    def test_basic_read(self, tmp_path: Path):
        """Read a valid long-format orthogroups TSV."""
        tsv = tmp_path / "orthogroups.tsv"
        tsv.write_text(
            "orthogroup_id\tgene_id\tspecies_id\n"
            "OG0000000\tg1\tsp1\n"
            "OG0000000\tg2\tsp1\n"
            "OG0000000\tg3\tsp2\n"
        )
        records = read_orthogroups_long_format(tsv)
        assert len(records) == 3
        assert records[0] == {"orthogroup_id": "OG0000000", "gene_id": "g1", "species_id": "sp1"}

    def test_empty_file_header_only(self, tmp_path: Path):
        """File with only header returns empty list."""
        tsv = tmp_path / "orthogroups.tsv"
        tsv.write_text("orthogroup_id\tgene_id\tspecies_id\n")
        records = read_orthogroups_long_format(tsv)
        assert records == []

    def test_missing_file(self, tmp_path: Path):
        """Missing file raises CompGeneError with E_INPUT_MISSING."""
        with pytest.raises(CompGeneError) as exc_info:
            read_orthogroups_long_format(tmp_path / "nonexistent.tsv")
        assert exc_info.value.error_code.value == "E_INPUT_MISSING"

    def test_bad_header(self, tmp_path: Path):
        """Invalid header raises CompGeneError with E_INPUT_FORMAT."""
        tsv = tmp_path / "bad.tsv"
        tsv.write_text("wrong_col1\twrong_col2\n")
        with pytest.raises(CompGeneError) as exc_info:
            read_orthogroups_long_format(tsv)
        assert exc_info.value.error_code.value == "E_INPUT_FORMAT"

    def test_strips_whitespace(self, tmp_path: Path):
        """Values with whitespace are stripped."""
        tsv = tmp_path / "orthogroups.tsv"
        tsv.write_text(
            "orthogroup_id\tgene_id\tspecies_id\n"
            "  OG0 \t g1 \t sp1 \n"
        )
        records = read_orthogroups_long_format(tsv)
        assert len(records) == 1
        assert records[0] == {"orthogroup_id": "OG0", "gene_id": "g1", "species_id": "sp1"}

    def test_skips_empty_fields(self, tmp_path: Path):
        """Rows with empty fields are skipped."""
        tsv = tmp_path / "orthogroups.tsv"
        tsv.write_text(
            "orthogroup_id\tgene_id\tspecies_id\n"
            "OG0\tg1\tsp1\n"
            "OG0\t\tsp2\n"
            "OG0\tg3\tsp2\n"
        )
        records = read_orthogroups_long_format(tsv)
        assert len(records) == 2


# =============================================================================
# Tests: End-to-end Integration (file → read → matrix → write → verify)
# =============================================================================

class TestEndToEndIntegration:
    def test_full_pipeline(self, tmp_path: Path):
        """Read long-format TSV, build both matrices, write and verify output."""
        # Write input file
        input_tsv = tmp_path / "orthogroups.tsv"
        input_tsv.write_text(
            "orthogroup_id\tgene_id\tspecies_id\n"
            "OG0000000\tg1\tspecies1\n"
            "OG0000000\tg2\tspecies1\n"
            "OG0000000\tg3\tspecies2\n"
            "OG0000001\tg4\tspecies1\n"
            "OG0000001\tg5\tspecies2\n"
            "OG0000002\tg6\tspecies2\n"
        )

        # Read
        records = read_orthogroups_long_format(input_tsv)
        assert len(records) == 6

        # Build matrices
        species = get_sorted_species(records)
        assert species == ["species1", "species2"]

        cn_matrix = build_copy_number_matrix(records)
        pa_matrix = build_presence_absence_matrix(records)

        # Verify CN values
        assert cn_matrix["OG0000000"]["species1"] == 2  # multi-copy
        assert cn_matrix["OG0000000"]["species2"] == 1
        assert cn_matrix["OG0000002"]["species1"] == 0  # absent
        assert cn_matrix["OG0000002"]["species2"] == 1

        # Verify PA values
        assert pa_matrix["OG0000000"]["species1"] == 1  # binary, not 2
        assert pa_matrix["OG0000002"]["species1"] == 0

        # Write and verify PA output
        pa_out = tmp_path / "pa.tsv"
        write_matrix_tsv(pa_matrix, species, pa_out)
        pa_lines = pa_out.read_text().strip().split("\n")
        assert pa_lines[0] == "orthogroup_id\tspecies1\tspecies2"
        assert pa_lines[1] == "OG0000000\t1\t1"
        assert pa_lines[2] == "OG0000001\t1\t1"
        assert pa_lines[3] == "OG0000002\t0\t1"

        # Write and verify CN output
        cn_out = tmp_path / "cn.tsv"
        write_matrix_tsv(cn_matrix, species, cn_out)
        cn_lines = cn_out.read_text().strip().split("\n")
        assert cn_lines[1] == "OG0000000\t2\t1"  # species1 has 2 copies
