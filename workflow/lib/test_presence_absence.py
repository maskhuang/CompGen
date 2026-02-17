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
    classify_orthogroup_sharing,
    get_sorted_species,
    read_orthogroups_long_format,
    read_presence_absence_matrix,
    summarize_sharing_counts,
    write_comparison_tsv,
    write_matrix_tsv,
    write_summary_tsv,
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


# =============================================================================
# Tests: read_presence_absence_matrix (Story 6A.2)
# =============================================================================

class TestReadPresenceAbsenceMatrix:
    def test_basic_read(self, tmp_path: Path):
        """Read a valid PA matrix TSV."""
        tsv = tmp_path / "pa.tsv"
        tsv.write_text(
            "orthogroup_id\tsp1\tsp2\tsp3\n"
            "OG0000000\t1\t1\t1\n"
            "OG0000001\t1\t0\t0\n"
            "OG0000002\t0\t1\t1\n"
        )
        matrix, species = read_presence_absence_matrix(tsv)
        assert species == ["sp1", "sp2", "sp3"]
        assert len(matrix) == 3
        assert matrix["OG0000000"] == {"sp1": 1, "sp2": 1, "sp3": 1}
        assert matrix["OG0000001"] == {"sp1": 1, "sp2": 0, "sp3": 0}
        assert matrix["OG0000002"] == {"sp1": 0, "sp2": 1, "sp3": 1}

    def test_missing_file(self, tmp_path: Path):
        """Missing file raises CompGeneError with E_INPUT_MISSING."""
        with pytest.raises(CompGeneError) as exc_info:
            read_presence_absence_matrix(tmp_path / "nonexistent.tsv")
        assert exc_info.value.error_code.value == "E_INPUT_MISSING"

    def test_bad_header(self, tmp_path: Path):
        """Header missing orthogroup_id raises CompGeneError."""
        tsv = tmp_path / "bad.tsv"
        tsv.write_text("wrong_col\tsp1\n1\t1\n")
        with pytest.raises(CompGeneError) as exc_info:
            read_presence_absence_matrix(tsv)
        assert exc_info.value.error_code.value == "E_INPUT_FORMAT"

    def test_empty_matrix(self, tmp_path: Path):
        """File with only header returns empty matrix."""
        tsv = tmp_path / "pa.tsv"
        tsv.write_text("orthogroup_id\tsp1\tsp2\n")
        matrix, species = read_presence_absence_matrix(tsv)
        assert matrix == {}
        assert species == ["sp1", "sp2"]

    def test_single_species(self, tmp_path: Path):
        """Single-species matrix read correctly."""
        tsv = tmp_path / "pa.tsv"
        tsv.write_text(
            "orthogroup_id\tspeciesA\n"
            "OG0000000\t1\n"
            "OG0000001\t1\n"
        )
        matrix, species = read_presence_absence_matrix(tsv)
        assert species == ["speciesA"]
        assert matrix["OG0000000"] == {"speciesA": 1}


# =============================================================================
# Tests: classify_orthogroup_sharing (Story 6A.2)
# =============================================================================

class TestClassifyOrthogroupSharing:
    def test_all_three_categories(self):
        """Matrix with all_shared, species_specific, and partial OGs."""
        matrix = {
            "OG0000000": {"sp1": 1, "sp2": 1, "sp3": 1},  # all_shared
            "OG0000001": {"sp1": 1, "sp2": 1, "sp3": 0},  # partial
            "OG0000002": {"sp1": 0, "sp2": 1, "sp3": 0},  # species_specific
        }
        species = ["sp1", "sp2", "sp3"]
        result = classify_orthogroup_sharing(matrix, species)

        assert len(result) == 3
        by_og = {r["orthogroup_id"]: r for r in result}

        assert by_og["OG0000000"]["category"] == "all_shared"
        assert by_og["OG0000000"]["n_species_present"] == 3
        assert by_og["OG0000000"]["present_species"] == "sp1,sp2,sp3"
        assert by_og["OG0000000"]["absent_species"] == ""
        assert by_og["OG0000000"]["specific_to"] == ""

        assert by_og["OG0000001"]["category"] == "partial"
        assert by_og["OG0000001"]["n_species_present"] == 2
        assert by_og["OG0000001"]["present_species"] == "sp1,sp2"
        assert by_og["OG0000001"]["absent_species"] == "sp3"

        assert by_og["OG0000002"]["category"] == "species_specific"
        assert by_og["OG0000002"]["n_species_present"] == 1
        assert by_og["OG0000002"]["specific_to"] == "sp2"

    def test_empty_matrix(self):
        """Empty matrix returns empty classifications."""
        result = classify_orthogroup_sharing({}, ["sp1", "sp2"])
        assert result == []

    def test_all_shared_only(self):
        """All orthogroups shared by all species."""
        matrix = {
            "OG0": {"sp1": 1, "sp2": 1},
            "OG1": {"sp1": 1, "sp2": 1},
        }
        result = classify_orthogroup_sharing(matrix, ["sp1", "sp2"])
        assert all(r["category"] == "all_shared" for r in result)

    def test_single_species_all_specific(self):
        """Single-species matrix: all OGs are species_specific."""
        matrix = {
            "OG0": {"sp1": 1},
            "OG1": {"sp1": 1},
        }
        result = classify_orthogroup_sharing(matrix, ["sp1"])
        assert all(r["category"] == "species_specific" for r in result)
        assert all(r["specific_to"] == "sp1" for r in result)

    def test_sorted_by_orthogroup_id(self):
        """Results are sorted by orthogroup_id."""
        matrix = {
            "OG0000002": {"sp1": 1, "sp2": 0},
            "OG0000000": {"sp1": 1, "sp2": 1},
            "OG0000001": {"sp1": 0, "sp2": 1},
        }
        result = classify_orthogroup_sharing(matrix, ["sp1", "sp2"])
        og_ids = [r["orthogroup_id"] for r in result]
        assert og_ids == ["OG0000000", "OG0000001", "OG0000002"]

    def test_two_species_no_partial(self):
        """With 2 species, partial requires exactly 2 present (== all_shared)."""
        matrix = {
            "OG0": {"sp1": 1, "sp2": 1},  # all_shared (2/2)
            "OG1": {"sp1": 1, "sp2": 0},  # species_specific (1/2)
            "OG2": {"sp1": 0, "sp2": 1},  # species_specific (1/2)
        }
        result = classify_orthogroup_sharing(matrix, ["sp1", "sp2"])
        by_og = {r["orthogroup_id"]: r for r in result}
        assert by_og["OG0"]["category"] == "all_shared"
        assert by_og["OG1"]["category"] == "species_specific"
        assert by_og["OG2"]["category"] == "species_specific"


# =============================================================================
# Tests: summarize_sharing_counts (Story 6A.2)
# =============================================================================

class TestSummarizeSharingCounts:
    def test_basic_summary(self):
        """Summary counts match classification categories."""
        classifications = [
            {"orthogroup_id": "OG0", "category": "all_shared", "specific_to": ""},
            {"orthogroup_id": "OG1", "category": "partial", "specific_to": ""},
            {"orthogroup_id": "OG2", "category": "species_specific", "specific_to": "sp1"},
            {"orthogroup_id": "OG3", "category": "species_specific", "specific_to": "sp2"},
            {"orthogroup_id": "OG4", "category": "species_specific", "specific_to": "sp1"},
        ]
        species = ["sp1", "sp2", "sp3"]
        result = summarize_sharing_counts(classifications, species)

        by_metric = {r["metric"]: r["count"] for r in result}
        assert by_metric["total_orthogroups"] == 5
        assert by_metric["all_shared"] == 1
        assert by_metric["partial_shared"] == 1
        assert by_metric["species_specific_total"] == 3
        assert by_metric["species_specific_sp1"] == 2
        assert by_metric["species_specific_sp2"] == 1
        assert by_metric["species_specific_sp3"] == 0

    def test_empty_classifications(self):
        """Empty classifications produce zero counts."""
        result = summarize_sharing_counts([], ["sp1", "sp2"])
        by_metric = {r["metric"]: r["count"] for r in result}
        assert by_metric["total_orthogroups"] == 0
        assert by_metric["all_shared"] == 0
        assert by_metric["partial_shared"] == 0
        assert by_metric["species_specific_total"] == 0

    def test_per_species_order(self):
        """Per-species rows follow species_list order."""
        result = summarize_sharing_counts([], ["zebra", "alpha", "mouse"])
        per_species_metrics = [
            r["metric"] for r in result
            if r["metric"].startswith("species_specific_") and r["metric"] != "species_specific_total"
        ]
        assert per_species_metrics == [
            "species_specific_zebra",
            "species_specific_alpha",
            "species_specific_mouse",
        ]


# =============================================================================
# Tests: write_comparison_tsv / write_summary_tsv (Story 6A.2)
# =============================================================================

class TestWriteComparisonTsv:
    def test_basic_write(self, tmp_path: Path):
        """Write classification results to TSV."""
        classifications = [
            {
                "orthogroup_id": "OG0",
                "category": "all_shared",
                "n_species_present": 3,
                "present_species": "sp1,sp2,sp3",
                "absent_species": "",
                "specific_to": "",
            },
            {
                "orthogroup_id": "OG1",
                "category": "species_specific",
                "n_species_present": 1,
                "present_species": "sp2",
                "absent_species": "sp1,sp3",
                "specific_to": "sp2",
            },
        ]
        output = tmp_path / "comparison.tsv"
        write_comparison_tsv(classifications, output)
        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 rows
        header = lines[0].split("\t")
        assert header == [
            "orthogroup_id", "category", "n_species_present",
            "present_species", "absent_species", "specific_to",
        ]
        row0 = lines[1].split("\t")
        assert row0[0] == "OG0"
        assert row0[1] == "all_shared"

    def test_atomic_write(self, tmp_path: Path):
        """No .tmp file remains after write."""
        output = tmp_path / "comparison.tsv"
        write_comparison_tsv([], output)
        assert not output.with_suffix(".tsv.tmp").exists()
        assert output.exists()

    def test_creates_parent_dirs(self, tmp_path: Path):
        """Parent directories are created."""
        output = tmp_path / "sub" / "dir" / "comparison.tsv"
        write_comparison_tsv([], output)
        assert output.exists()


class TestWriteSummaryTsv:
    def test_basic_write(self, tmp_path: Path):
        """Write summary rows to TSV."""
        rows = [
            {"metric": "total_orthogroups", "count": 5},
            {"metric": "all_shared", "count": 2},
        ]
        output = tmp_path / "summary.tsv"
        write_summary_tsv(rows, output)
        assert output.exists()
        lines = output.read_text().strip().split("\n")
        assert len(lines) == 3
        assert lines[0] == "metric\tcount"
        assert lines[1] == "total_orthogroups\t5"
        assert lines[2] == "all_shared\t2"

    def test_atomic_write(self, tmp_path: Path):
        """No .tmp file remains after write."""
        output = tmp_path / "summary.tsv"
        write_summary_tsv([], output)
        assert not output.with_suffix(".tsv.tmp").exists()


# =============================================================================
# Tests: End-to-end Story 6A.2 Integration
# =============================================================================

class TestComparisonEndToEnd:
    def test_full_pipeline(self, tmp_path: Path):
        """Read PA matrix → classify → summarize → write → verify."""
        # Write input PA matrix
        pa_tsv = tmp_path / "presence_absence.tsv"
        pa_tsv.write_text(
            "orthogroup_id\tspecies1\tspecies2\tspecies3\n"
            "OG0000000\t1\t1\t1\n"   # all_shared
            "OG0000001\t1\t1\t0\n"   # partial
            "OG0000002\t0\t1\t0\n"   # species_specific (species2)
            "OG0000003\t1\t0\t0\n"   # species_specific (species1)
            "OG0000004\t0\t1\t1\n"   # partial
        )

        # Read PA matrix
        matrix, species = read_presence_absence_matrix(pa_tsv)
        assert len(matrix) == 5
        assert species == ["species1", "species2", "species3"]

        # Classify
        classifications = classify_orthogroup_sharing(matrix, species)
        assert len(classifications) == 5
        by_og = {c["orthogroup_id"]: c for c in classifications}
        assert by_og["OG0000000"]["category"] == "all_shared"
        assert by_og["OG0000001"]["category"] == "partial"
        assert by_og["OG0000002"]["category"] == "species_specific"
        assert by_og["OG0000002"]["specific_to"] == "species2"
        assert by_og["OG0000003"]["category"] == "species_specific"
        assert by_og["OG0000003"]["specific_to"] == "species1"
        assert by_og["OG0000004"]["category"] == "partial"

        # Summarize
        summary = summarize_sharing_counts(classifications, species)
        by_metric = {r["metric"]: r["count"] for r in summary}
        assert by_metric["total_orthogroups"] == 5
        assert by_metric["all_shared"] == 1
        assert by_metric["partial_shared"] == 2
        assert by_metric["species_specific_total"] == 2
        assert by_metric["species_specific_species1"] == 1
        assert by_metric["species_specific_species2"] == 1
        assert by_metric["species_specific_species3"] == 0

        # Write and verify comparison
        comp_out = tmp_path / "comparison.tsv"
        write_comparison_tsv(classifications, comp_out)
        comp_lines = comp_out.read_text().strip().split("\n")
        assert len(comp_lines) == 6  # header + 5 OGs
        assert comp_lines[0].startswith("orthogroup_id\t")

        # Write and verify summary
        sum_out = tmp_path / "summary.tsv"
        write_summary_tsv(summary, sum_out)
        sum_lines = sum_out.read_text().strip().split("\n")
        assert sum_lines[0] == "metric\tcount"
        assert "total_orthogroups\t5" in sum_lines[1]
