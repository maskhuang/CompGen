"""
Unit tests for LiftoffAdapter.

Source: Story 5.1 Liftoff Adapter
"""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import subprocess
import tempfile

from workflow.adapters.liftoff import (
    LiftoffAdapter,
    parse_liftoff_stats,
    parse_liftoff_gff_coverage,
)
from workflow.adapters.base import AdapterContext, ToolSpec
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def adapter():
    """Create a LiftoffAdapter instance."""
    return LiftoffAdapter()


@pytest.fixture
def mock_context():
    """Create a mock AdapterContext for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create mock input files
        ref_gff = tmpdir / "reference.gff3"
        ref_fa = tmpdir / "reference.fa"
        target_fa = tmpdir / "target.fa"

        ref_gff.write_text("##gff-version 3\nchr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene1")
        ref_fa.write_text(">chr1\nATGC")
        target_fa.write_text(">chr1\nATGC")

        # Create output paths
        out_dir = tmpdir / "output"
        out_dir.mkdir()

        ctx = AdapterContext(
            inputs={
                "reference_gff": ref_gff,
                "reference_fa": ref_fa,
                "target_fa": target_fa,
            },
            outputs={
                "lifted_gff": out_dir / "lifted.gff3",
                "unmapped": out_dir / "unmapped.txt",
                "intermediate_dir": out_dir / "intermediate",
            },
            config={
                "liftoff": {
                    "min_coverage": 0.5,
                    "min_identity": 0.5,
                    "timeout": 3600,
                }
            },
            wildcards={"reference": "human", "target": "mouse"},
            threads=4,
        )

        yield ctx


@pytest.fixture
def sample_gff_content():
    """Sample GFF content with coverage attributes."""
    return """##gff-version 3
chr1\tliftoff\tgene\t1\t1000\t.\t+\t.\tID=gene1;coverage=0.95
chr1\tliftoff\tmRNA\t1\t1000\t.\t+\t.\tID=mrna1;Parent=gene1;coverage=0.95
chr1\tliftoff\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=mrna1;coverage=0.98
chr1\tliftoff\texon\t600\t1000\t.\t+\t.\tID=exon2;Parent=mrna1;coverage=0.90
chr1\tliftoff\tCDS\t100\t900\t.\t+\t0\tID=cds1;Parent=mrna1;coverage=0.92
"""


# =============================================================================
# ToolSpec Tests
# =============================================================================

class TestLiftoffSpec:
    """Tests for LiftoffAdapter.spec property."""

    def test_spec_name(self, adapter):
        """Test that spec name is 'liftoff'."""
        assert adapter.spec.name == "liftoff"

    def test_spec_min_version(self, adapter):
        """Test that spec has minimum version 1.6.0."""
        assert adapter.spec.min_version == "1.6.0"

    def test_spec_max_version(self, adapter):
        """Test that spec has maximum version constraint."""
        assert adapter.spec.max_version == "1.99.99"

    def test_spec_conda_env(self, adapter):
        """Test that spec references liftoff.yaml conda environment."""
        assert adapter.spec.conda_env == "liftoff.yaml"

    def test_spec_version_compatible_1_6_3(self, adapter):
        """Test that version 1.6.3 is compatible."""
        assert adapter.spec.check_version_compatible("1.6.3") is True

    def test_spec_version_compatible_1_6_0(self, adapter):
        """Test that version 1.6.0 is compatible."""
        assert adapter.spec.check_version_compatible("1.6.0") is True

    def test_spec_version_incompatible_1_5_0(self, adapter):
        """Test that version 1.5.0 is incompatible."""
        assert adapter.spec.check_version_compatible("1.5.0") is False

    def test_spec_version_incompatible_2_0_0(self, adapter):
        """Test that version 2.0.0 is incompatible."""
        assert adapter.spec.check_version_compatible("2.0.0") is False


# =============================================================================
# Version Check Tests
# =============================================================================

class TestCheckVersion:
    """Tests for LiftoffAdapter.check_version method."""

    def test_check_version_success(self, adapter):
        """Test successful version detection."""
        mock_result = MagicMock()
        mock_result.stdout = "liftoff 1.6.3\n"
        mock_result.stderr = ""

        with patch("subprocess.run", return_value=mock_result):
            version = adapter.check_version()
            assert version == "1.6.3"

    def test_check_version_stderr_output(self, adapter):
        """Test version detection from stderr."""
        mock_result = MagicMock()
        mock_result.stdout = ""
        mock_result.stderr = "liftoff 1.6.3\n"

        with patch("subprocess.run", return_value=mock_result):
            version = adapter.check_version()
            assert version == "1.6.3"

    def test_check_version_not_found(self, adapter):
        """Test error when Liftoff not installed."""
        with patch("subprocess.run", side_effect=FileNotFoundError):
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND

    def test_check_version_timeout(self, adapter):
        """Test error when version check times out."""
        with patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="liftoff", timeout=30)):
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND

    def test_check_version_incompatible(self, adapter):
        """Test error for unsupported version."""
        mock_result = MagicMock()
        mock_result.stdout = "liftoff 1.5.0\n"
        mock_result.stderr = ""

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION

    def test_check_version_parse_error(self, adapter):
        """Test error when version cannot be parsed."""
        mock_result = MagicMock()
        mock_result.stdout = "unknown output format"
        mock_result.stderr = ""

        with patch("subprocess.run", return_value=mock_result):
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND


# =============================================================================
# Input Validation Tests
# =============================================================================

class TestValidateInputs:
    """Tests for LiftoffAdapter.validate_inputs method."""

    def test_validate_inputs_success(self, adapter, mock_context):
        """Test successful input validation."""
        # Should not raise
        adapter.validate_inputs(mock_context)

    def test_validate_inputs_missing_ref_gff(self, adapter, mock_context):
        """Test error when reference GFF is missing."""
        mock_context.inputs["reference_gff"] = Path("/nonexistent/file.gff3")

        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(mock_context)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_validate_inputs_missing_ref_fa(self, adapter, mock_context):
        """Test error when reference FASTA is missing."""
        mock_context.inputs["reference_fa"] = Path("/nonexistent/file.fa")

        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(mock_context)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_validate_inputs_missing_target_fa(self, adapter, mock_context):
        """Test error when target FASTA is missing."""
        mock_context.inputs["target_fa"] = Path("/nonexistent/file.fa")

        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(mock_context)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING

    def test_validate_inputs_none_ref_gff(self, adapter, mock_context):
        """Test error when reference GFF is None."""
        del mock_context.inputs["reference_gff"]

        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(mock_context)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING


# =============================================================================
# Command Building Tests
# =============================================================================

class TestBuildCommand:
    """Tests for LiftoffAdapter.build_command method."""

    def test_build_command_basic(self, adapter, mock_context):
        """Test basic command construction."""
        cmd = adapter.build_command(mock_context)

        assert cmd[0] == "liftoff"
        assert "-g" in cmd
        assert "-o" in cmd
        assert "-u" in cmd
        assert "-dir" in cmd
        assert "-p" in cmd
        assert "4" in cmd  # threads

    def test_build_command_with_coverage(self, adapter, mock_context):
        """Test command with coverage threshold."""
        cmd = adapter.build_command(mock_context)

        assert "-a" in cmd
        assert "0.5" in cmd

    def test_build_command_with_identity(self, adapter, mock_context):
        """Test command with identity threshold."""
        cmd = adapter.build_command(mock_context)

        assert "-s" in cmd
        assert "0.5" in cmd

    def test_build_command_with_copies(self, adapter, mock_context):
        """Test command with copies option enabled."""
        mock_context.config["liftoff"]["copies"] = True
        cmd = adapter.build_command(mock_context)

        assert "-copies" in cmd

    def test_build_command_with_flank(self, adapter, mock_context):
        """Test command with flank option."""
        mock_context.config["liftoff"]["flank"] = 5000
        cmd = adapter.build_command(mock_context)

        assert "-flank" in cmd
        flank_idx = cmd.index("-flank")
        assert cmd[flank_idx + 1] == "5000"

    def test_build_command_flank_zero_not_added(self, adapter, mock_context):
        """Test that flank=0 does not add -flank flag."""
        mock_context.config["liftoff"]["flank"] = 0
        cmd = adapter.build_command(mock_context)

        assert "-flank" not in cmd

    def test_build_command_positional_args(self, adapter, mock_context):
        """Test that target and reference FASTAs are at end."""
        cmd = adapter.build_command(mock_context)

        # Last two arguments should be target and reference FASTAs
        assert "target.fa" in str(cmd[-2])
        assert "reference.fa" in str(cmd[-1])


# =============================================================================
# Output Parsing Tests
# =============================================================================

class TestParseOutputs:
    """Tests for LiftoffAdapter.parse_outputs method."""

    def test_parse_outputs_success(self, adapter, mock_context, sample_gff_content):
        """Test successful output parsing."""
        # Create output files
        mock_context.outputs["lifted_gff"].write_text(sample_gff_content)
        mock_context.outputs["unmapped"].write_text("gene_x\ngene_y\n")

        result = adapter.parse_outputs(mock_context)

        assert "lifted_gff" in result.outputs
        assert "unmapped" in result.outputs
        assert result.summary["lifted_genes"] == 1  # Only 1 gene in sample
        assert result.summary["lifted_features"] == 5  # All feature lines
        assert result.summary["unmapped_genes"] == 2
        assert result.summary["lift_rate"] > 0

    def test_parse_outputs_empty(self, adapter, mock_context):
        """Test output parsing with empty files."""
        mock_context.outputs["lifted_gff"].write_text("")
        mock_context.outputs["unmapped"].write_text("")

        result = adapter.parse_outputs(mock_context)

        assert result.summary["lifted_genes"] == 0
        assert result.summary["lifted_features"] == 0
        assert result.summary["unmapped_genes"] == 0
        assert result.summary["lift_rate"] == 0.0


# =============================================================================
# Timeout Tests
# =============================================================================

class TestTimeout:
    """Tests for LiftoffAdapter.timeout_seconds method."""

    def test_timeout_default(self, adapter, mock_context):
        """Test default timeout of 60 minutes."""
        mock_context.config["liftoff"].pop("timeout", None)
        timeout = adapter.timeout_seconds(mock_context)
        assert timeout == 3600

    def test_timeout_configured(self, adapter, mock_context):
        """Test configured timeout."""
        mock_context.config["liftoff"]["timeout"] = 7200
        timeout = adapter.timeout_seconds(mock_context)
        assert timeout == 7200


# =============================================================================
# Error Classification Tests
# =============================================================================

class TestClassifyError:
    """Tests for LiftoffAdapter.classify_error method."""

    def test_classify_file_not_found(self, adapter, mock_context):
        """Test classification of file not found error."""
        error_code, retryable = adapter.classify_error(
            mock_context, 1, "Error: No such file or directory"
        )
        assert error_code == ErrorCode.E_INPUT_MISSING
        assert retryable is False

    def test_classify_gff_format_error(self, adapter, mock_context):
        """Test classification of GFF format error."""
        error_code, retryable = adapter.classify_error(
            mock_context, 1, "Invalid GFF format error"
        )
        assert error_code == ErrorCode.E_INPUT_FORMAT
        assert retryable is False

    def test_classify_fasta_format_error(self, adapter, mock_context):
        """Test classification of FASTA format error."""
        error_code, retryable = adapter.classify_error(
            mock_context, 1, "FASTA parse error"
        )
        assert error_code == ErrorCode.E_INPUT_FORMAT
        assert retryable is False

    def test_classify_memory_error(self, adapter, mock_context):
        """Test classification of memory error."""
        error_code, retryable = adapter.classify_error(
            mock_context, 1, "Out of memory"
        )
        assert error_code == ErrorCode.E_OOM
        assert retryable is False

    def test_classify_sigkill(self, adapter, mock_context):
        """Test classification of SIGKILL (timeout/OOM killer)."""
        error_code, retryable = adapter.classify_error(
            mock_context, -9, ""
        )
        assert error_code == ErrorCode.E_TIMEOUT
        assert retryable is True


# =============================================================================
# Helper Function Tests
# =============================================================================

class TestParseLiftoffStats:
    """Tests for parse_liftoff_stats function."""

    def test_parse_stats_basic(self, tmp_path, sample_gff_content):
        """Test basic stats parsing."""
        gff_path = tmp_path / "lifted.gff3"
        unmapped_path = tmp_path / "unmapped.txt"

        gff_path.write_text(sample_gff_content)
        unmapped_path.write_text("gene1\ngene2\ngene3\n")

        stats = parse_liftoff_stats(gff_path, unmapped_path)

        assert stats["lifted_genes"] == 1  # Only 1 gene in sample_gff_content
        assert stats["lifted_features"] == 5  # All feature lines
        assert stats["unmapped_genes"] == 3
        assert stats["total_genes"] == 4
        assert stats["lift_rate"] == pytest.approx(0.25, rel=0.01)  # 1/4

    def test_parse_stats_empty_files(self, tmp_path):
        """Test stats parsing with empty files."""
        gff_path = tmp_path / "lifted.gff3"
        unmapped_path = tmp_path / "unmapped.txt"

        gff_path.write_text("")
        unmapped_path.write_text("")

        stats = parse_liftoff_stats(gff_path, unmapped_path)

        assert stats["lifted_genes"] == 0
        assert stats["lifted_features"] == 0
        assert stats["unmapped_genes"] == 0
        assert stats["lift_rate"] == 0.0

    def test_parse_stats_missing_files(self, tmp_path):
        """Test stats parsing with missing files."""
        gff_path = tmp_path / "nonexistent.gff3"
        unmapped_path = tmp_path / "nonexistent.txt"

        stats = parse_liftoff_stats(gff_path, unmapped_path)

        assert stats["lifted_genes"] == 0
        assert stats["lifted_features"] == 0
        assert stats["unmapped_genes"] == 0


class TestParseLiftoffGffCoverage:
    """Tests for parse_liftoff_gff_coverage function."""

    def test_parse_coverage_basic(self, tmp_path, sample_gff_content):
        """Test basic coverage parsing."""
        gff_path = tmp_path / "lifted.gff3"
        gff_path.write_text(sample_gff_content)

        coverage = parse_liftoff_gff_coverage(gff_path)

        assert coverage["feature_count"] == 5
        assert coverage["min_coverage"] is not None
        assert coverage["max_coverage"] is not None
        assert coverage["mean_coverage"] is not None

    def test_parse_coverage_no_attributes(self, tmp_path):
        """Test coverage parsing without coverage attributes."""
        gff_path = tmp_path / "lifted.gff3"
        gff_path.write_text("chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene1")

        coverage = parse_liftoff_gff_coverage(gff_path)

        assert coverage["feature_count"] == 0
        assert coverage["mean_coverage"] == 0.0

    def test_parse_coverage_missing_file(self, tmp_path):
        """Test coverage parsing with missing file."""
        gff_path = tmp_path / "nonexistent.gff3"

        coverage = parse_liftoff_gff_coverage(gff_path)

        assert coverage["feature_count"] == 0
        assert coverage["mean_coverage"] == 0.0
