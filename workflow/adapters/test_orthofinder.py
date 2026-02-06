"""
Unit tests for workflow/adapters/orthofinder.py

Tests for OrthoFinderAdapter: OrthoFinder orthology inference tool integration.
"""

from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from workflow.adapters.base import (
    ToolSpec,
    RunResult,
    AdapterContext,
)
from workflow.lib.errors import ErrorCode, CompGeneError

from workflow.adapters.orthofinder import (
    OrthoFinderAdapter,
    parse_gene_count_tsv,
    count_unassigned_genes,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def adapter() -> OrthoFinderAdapter:
    """Create an OrthoFinderAdapter instance."""
    return OrthoFinderAdapter()


@pytest.fixture
def proteins_dir(tmp_path: Path) -> Path:
    """Create a temporary proteins directory with test .fa files."""
    pdir = tmp_path / "input_proteins"
    pdir.mkdir()

    (pdir / "species1.fa").write_text(">gene1\nMKTAYI\n>gene2\nMKTAYA\n")
    (pdir / "species2.fa").write_text(">gene3\nMKTAYI\n")
    (pdir / "species3.fa").write_text(">gene4\nMKTAYA\n>gene5\nMKTAYI\n>gene6\nMKT\n")

    return pdir


@pytest.fixture
def basic_context(tmp_path: Path, proteins_dir: Path) -> AdapterContext:
    """Create a basic AdapterContext for testing."""
    results_dir = tmp_path / "results" / "orthology"
    results_dir.mkdir(parents=True)

    return AdapterContext(
        inputs={"proteins_dir": proteins_dir},
        outputs={
            "results_dir": results_dir,
            "output_dir": results_dir,
        },
        config={
            "orthofinder": {
                "search_method": "diamond",
                "timeout": 14400,
            }
        },
        wildcards={},
        threads=8,
    )


@pytest.fixture
def gene_count_file(tmp_path: Path) -> Path:
    """Create a sample Orthogroups.GeneCount.tsv file."""
    gc_dir = tmp_path / "Orthogroups"
    gc_dir.mkdir(parents=True)
    gc_path = gc_dir / "Orthogroups.GeneCount.tsv"
    gc_path.write_text(
        "Orthogroup\tspecies1\tspecies2\tspecies3\tTotal\n"
        "OG0000000\t2\t1\t3\t6\n"
        "OG0000001\t1\t1\t1\t3\n"
        "OG0000002\t0\t1\t1\t2\n"
        "OG0000003\t1\t1\t1\t3\n"
    )
    return gc_path


@pytest.fixture
def unassigned_file(tmp_path: Path) -> Path:
    """Create a sample Orthogroups_UnassignedGenes.tsv file."""
    og_dir = tmp_path / "Orthogroups"
    og_dir.mkdir(parents=True, exist_ok=True)
    ua_path = og_dir / "Orthogroups_UnassignedGenes.tsv"
    ua_path.write_text(
        "Orthogroup\tspecies1\tspecies2\tspecies3\n"
        "OG0010000\tgene99\t\t\n"
        "OG0010001\t\tgene100\t\n"
    )
    return ua_path


# =============================================================================
# Tests: ToolSpec (Subtask 1.2)
# =============================================================================

class TestOrthoFinderSpec:
    """Tests for OrthoFinderAdapter.spec property."""

    def test_spec_name(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.name == "orthofinder"

    def test_spec_min_version(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.min_version == "2.5.0"

    def test_spec_max_version(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.max_version == "2.5.99"

    def test_spec_conda_env(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.conda_env == "orthofinder.yaml"

    def test_spec_version_compatible_255(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.check_version_compatible("2.5.5") is True

    def test_spec_version_compatible_250(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.check_version_compatible("2.5.0") is True

    def test_spec_version_incompatible_240(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.check_version_compatible("2.4.0") is False

    def test_spec_version_incompatible_300(self, adapter: OrthoFinderAdapter):
        assert adapter.spec.check_version_compatible("3.0.0") is False


# =============================================================================
# Tests: check_version (Subtask 1.3)
# =============================================================================

class TestCheckVersion:
    """Tests for OrthoFinderAdapter.check_version()."""

    @patch("subprocess.run")
    def test_version_255(self, mock_run, adapter: OrthoFinderAdapter):
        mock_run.return_value = MagicMock(
            stdout="OrthoFinder version 2.5.5\n",
            stderr="",
        )
        assert adapter.check_version() == "2.5.5"

    @patch("subprocess.run")
    def test_version_250(self, mock_run, adapter: OrthoFinderAdapter):
        mock_run.return_value = MagicMock(
            stdout="OrthoFinder version 2.5.0\n",
            stderr="",
        )
        assert adapter.check_version() == "2.5.0"

    @patch("subprocess.run")
    def test_version_in_stderr(self, mock_run, adapter: OrthoFinderAdapter):
        mock_run.return_value = MagicMock(
            stdout="",
            stderr="OrthoFinder version 2.5.5\n",
        )
        assert adapter.check_version() == "2.5.5"

    @patch("subprocess.run")
    def test_version_unsupported_240(self, mock_run, adapter: OrthoFinderAdapter):
        mock_run.return_value = MagicMock(
            stdout="OrthoFinder version 2.4.0\n",
            stderr="",
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.check_version()
        assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION

    @patch("subprocess.run")
    def test_version_unsupported_300(self, mock_run, adapter: OrthoFinderAdapter):
        mock_run.return_value = MagicMock(
            stdout="OrthoFinder version 3.0.0\n",
            stderr="",
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.check_version()
        assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION

    @patch("subprocess.run")
    def test_version_unparseable(self, mock_run, adapter: OrthoFinderAdapter):
        mock_run.return_value = MagicMock(
            stdout="unknown output\n",
            stderr="",
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.check_version()
        assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND

    @patch("subprocess.run", side_effect=FileNotFoundError())
    def test_version_not_installed(self, mock_run, adapter: OrthoFinderAdapter):
        with pytest.raises(CompGeneError) as exc_info:
            adapter.check_version()
        assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND

    @patch("subprocess.run", side_effect=TimeoutError())
    def test_version_timeout(self, mock_run, adapter: OrthoFinderAdapter):
        """TimeoutExpired inherits from SubprocessError, not TimeoutError."""
        import subprocess as sp
        with patch("subprocess.run", side_effect=sp.TimeoutExpired(cmd="orthofinder", timeout=30)):
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND


# =============================================================================
# Tests: validate_inputs (Subtask 1.4)
# =============================================================================

class TestValidateInputs:
    """Tests for OrthoFinderAdapter.validate_inputs()."""

    def test_valid_inputs(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        # Should not raise
        adapter.validate_inputs(basic_context)

    def test_missing_proteins_dir_key(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},
            wildcards={},
            threads=1,
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_nonexistent_proteins_dir(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        ctx = AdapterContext(
            inputs={"proteins_dir": tmp_path / "nonexistent"},
            outputs={},
            config={},
            wildcards={},
            threads=1,
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_insufficient_species(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        pdir = tmp_path / "proteins"
        pdir.mkdir()
        (pdir / "only_one.fa").write_text(">gene1\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins_dir": pdir},
            outputs={},
            config={},
            wildcards={},
            threads=1,
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING
        assert "at least 2" in str(exc_info.value.message)

    def test_empty_directory(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        pdir = tmp_path / "empty"
        pdir.mkdir()

        ctx = AdapterContext(
            inputs={"proteins_dir": pdir},
            outputs={},
            config={},
            wildcards={},
            threads=1,
        )
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_exactly_two_species(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        pdir = tmp_path / "two_species"
        pdir.mkdir()
        (pdir / "sp1.fa").write_text(">g1\nMKT\n")
        (pdir / "sp2.fa").write_text(">g2\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins_dir": pdir},
            outputs={},
            config={},
            wildcards={},
            threads=1,
        )
        # Should not raise
        adapter.validate_inputs(ctx)


# =============================================================================
# Tests: build_command (Subtask 1.5)
# =============================================================================

class TestBuildCommand:
    """Tests for OrthoFinderAdapter.build_command()."""

    def test_basic_command(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        cmd = adapter.build_command(basic_context)
        assert cmd[0] == "orthofinder"
        assert "-f" in cmd
        assert "-t" in cmd
        assert "-a" in cmd
        assert "-S" in cmd
        assert "-n" in cmd

    def test_threads_parameter(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        cmd = adapter.build_command(basic_context)
        t_idx = cmd.index("-t")
        assert cmd[t_idx + 1] == "8"

    def test_analysis_threads_capped(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        """Analysis threads should be capped at 4."""
        cmd = adapter.build_command(basic_context)
        a_idx = cmd.index("-a")
        assert int(cmd[a_idx + 1]) <= 4

    def test_search_method_diamond(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        cmd = adapter.build_command(basic_context)
        s_idx = cmd.index("-S")
        assert cmd[s_idx + 1] == "diamond"

    def test_search_method_blast(self, adapter: OrthoFinderAdapter, proteins_dir: Path, tmp_path: Path):
        ctx = AdapterContext(
            inputs={"proteins_dir": proteins_dir},
            outputs={"output_dir": tmp_path / "output"},
            config={"orthofinder": {"search_method": "blast"}},
            wildcards={},
            threads=4,
        )
        cmd = adapter.build_command(ctx)
        s_idx = cmd.index("-S")
        assert cmd[s_idx + 1] == "blast"

    def test_output_dir_included(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        cmd = adapter.build_command(basic_context)
        assert "-o" in cmd

    def test_result_suffix(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        cmd = adapter.build_command(basic_context)
        n_idx = cmd.index("-n")
        assert cmd[n_idx + 1] == "compgene"

    def test_default_search_method(self, adapter: OrthoFinderAdapter, proteins_dir: Path, tmp_path: Path):
        """Default search method should be diamond when not configured."""
        ctx = AdapterContext(
            inputs={"proteins_dir": proteins_dir},
            outputs={"output_dir": tmp_path / "output"},
            config={},
            wildcards={},
            threads=4,
        )
        cmd = adapter.build_command(ctx)
        s_idx = cmd.index("-S")
        assert cmd[s_idx + 1] == "diamond"

    def test_low_thread_count(self, adapter: OrthoFinderAdapter, proteins_dir: Path, tmp_path: Path):
        """Analysis threads should equal search threads when threads < 4."""
        ctx = AdapterContext(
            inputs={"proteins_dir": proteins_dir},
            outputs={"output_dir": tmp_path / "output"},
            config={},
            wildcards={},
            threads=2,
        )
        cmd = adapter.build_command(ctx)
        a_idx = cmd.index("-a")
        assert cmd[a_idx + 1] == "2"


# =============================================================================
# Tests: expected_outputs (Subtask 1.6)
# =============================================================================

class TestExpectedOutputs:
    """Tests for OrthoFinderAdapter.expected_outputs()."""

    def test_expected_files(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        outputs = adapter.expected_outputs(basic_context)
        assert len(outputs) == 3

        output_names = [p.name for p in outputs]
        assert "Orthogroups.tsv" in output_names
        assert "Orthogroups.GeneCount.tsv" in output_names
        assert "SpeciesTree_rooted.txt" in output_names

    def test_output_paths_under_results_dir(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        results_dir = basic_context.outputs["results_dir"]
        outputs = adapter.expected_outputs(basic_context)
        for output in outputs:
            assert str(output).startswith(str(results_dir))


# =============================================================================
# Tests: parse_outputs (Subtask 1.7)
# =============================================================================

class TestParseOutputs:
    """Tests for OrthoFinderAdapter.parse_outputs()."""

    def test_parse_with_files(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        results_dir = tmp_path / "results"
        og_dir = results_dir / "Orthogroups"
        st_dir = results_dir / "Species_Tree"
        og_dir.mkdir(parents=True)
        st_dir.mkdir(parents=True)

        # Create test output files
        (og_dir / "Orthogroups.tsv").write_text(
            "Orthogroup\tsp1\tsp2\n"
            "OG0000000\tg1\tg2\n"
        )
        (og_dir / "Orthogroups.GeneCount.tsv").write_text(
            "Orthogroup\tsp1\tsp2\tTotal\n"
            "OG0000000\t1\t1\t2\n"
            "OG0000001\t2\t1\t3\n"
            "OG0000002\t1\t1\t2\n"
        )
        (og_dir / "Orthogroups_UnassignedGenes.tsv").write_text(
            "Orthogroup\tsp1\tsp2\n"
            "OG0010000\tgene99\t\n"
        )
        (st_dir / "SpeciesTree_rooted.txt").write_text("(sp1,sp2);")

        ctx = AdapterContext(
            inputs={},
            outputs={"results_dir": results_dir},
            config={},
            wildcards={},
            threads=1,
        )

        result = adapter.parse_outputs(ctx)
        assert isinstance(result, RunResult)
        assert result.summary["total_orthogroups"] == 3
        assert result.summary["species_count"] == 2
        assert result.summary["single_copy_orthogroups"] == 2
        assert result.summary["unassigned_genes"] == 1
        assert "sp1" in result.summary["species_names"]
        assert "sp2" in result.summary["species_names"]

    def test_parse_empty_outputs(self, adapter: OrthoFinderAdapter, tmp_path: Path):
        results_dir = tmp_path / "results"
        results_dir.mkdir(parents=True)

        ctx = AdapterContext(
            inputs={},
            outputs={"results_dir": results_dir},
            config={},
            wildcards={},
            threads=1,
        )

        result = adapter.parse_outputs(ctx)
        assert result.summary["total_orthogroups"] == 0
        assert result.summary["unassigned_genes"] == 0


# =============================================================================
# Tests: timeout_seconds (Subtask 1.8)
# =============================================================================

class TestTimeoutSeconds:
    """Tests for OrthoFinderAdapter.timeout_seconds()."""

    def test_default_timeout(self, adapter: OrthoFinderAdapter):
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},
            wildcards={},
            threads=1,
        )
        assert adapter.timeout_seconds(ctx) == 14400  # 4 hours

    def test_configured_timeout(self, adapter: OrthoFinderAdapter):
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={"orthofinder": {"timeout": 7200}},
            wildcards={},
            threads=1,
        )
        assert adapter.timeout_seconds(ctx) == 7200

    def test_missing_orthofinder_config(self, adapter: OrthoFinderAdapter):
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={"other": {}},
            wildcards={},
            threads=1,
        )
        assert adapter.timeout_seconds(ctx) == 14400


# =============================================================================
# Tests: classify_error (Subtask 1.9)
# =============================================================================

class TestClassifyError:
    """Tests for OrthoFinderAdapter.classify_error()."""

    def test_diamond_error(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 1, "DIAMOND error: database not found"
        )
        assert code == ErrorCode.E_NONZERO_EXIT
        assert retryable is False

    def test_no_fasta_files(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 1, "ERROR: No fasta files found in directory"
        )
        assert code == ErrorCode.E_INPUT_MISSING
        assert retryable is False

    def test_no_sequences(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 1, "No sequences in input files"
        )
        assert code == ErrorCode.E_INPUT_MISSING
        assert retryable is False

    def test_mcl_error(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 1, "MCL error during clustering"
        )
        assert code == ErrorCode.E_NONZERO_EXIT
        assert retryable is False

    def test_memory_error(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 1, "Out of memory during execution"
        )
        assert code == ErrorCode.E_OOM
        assert retryable is False

    def test_sigkill(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, -9, ""
        )
        assert code == ErrorCode.E_TIMEOUT
        assert retryable is True

    def test_generic_error(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 1, "some unknown error"
        )
        assert code == ErrorCode.E_NONZERO_EXIT
        assert retryable is False

    def test_command_not_found(self, adapter: OrthoFinderAdapter, basic_context: AdapterContext):
        code, retryable = adapter.classify_error(
            basic_context, 127, "command not found"
        )
        assert code == ErrorCode.E_TOOL_NOT_FOUND
        assert retryable is False


# =============================================================================
# Tests: Helper Functions
# =============================================================================

class TestParseGeneCountTsv:
    """Tests for parse_gene_count_tsv helper function."""

    def test_parse_valid_file(self, gene_count_file: Path):
        stats = parse_gene_count_tsv(gene_count_file)
        assert stats["total_orthogroups"] == 4
        assert stats["species_names"] == ["species1", "species2", "species3"]
        assert stats["single_copy_count"] == 2  # OG0000001 and OG0000003

    def test_parse_nonexistent_file(self, tmp_path: Path):
        stats = parse_gene_count_tsv(tmp_path / "nonexistent.tsv")
        assert stats["total_orthogroups"] == 0
        assert stats["species_names"] == []
        assert stats["single_copy_count"] == 0

    def test_parse_header_only(self, tmp_path: Path):
        f = tmp_path / "header_only.tsv"
        f.write_text("Orthogroup\tsp1\tsp2\tTotal\n")
        stats = parse_gene_count_tsv(f)
        assert stats["total_orthogroups"] == 0
        assert stats["species_names"] == ["sp1", "sp2"]


class TestCountUnassignedGenes:
    """Tests for count_unassigned_genes helper function."""

    def test_count_valid_file(self, unassigned_file: Path):
        count = count_unassigned_genes(unassigned_file)
        assert count == 2

    def test_count_nonexistent_file(self, tmp_path: Path):
        count = count_unassigned_genes(tmp_path / "nonexistent.tsv")
        assert count == 0

    def test_count_header_only(self, tmp_path: Path):
        f = tmp_path / "header.tsv"
        f.write_text("Orthogroup\tsp1\tsp2\n")
        count = count_unassigned_genes(f)
        assert count == 0
