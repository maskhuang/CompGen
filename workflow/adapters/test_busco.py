"""
Unit tests for workflow/adapters/busco.py

Tests for BuscoAdapter: BUSCO quality control tool integration.
"""

import re
from pathlib import Path
from typing import Any
from unittest.mock import Mock, patch, MagicMock

import pytest

from workflow.adapters.base import (
    ToolSpec,
    RunResult,
    AdapterContext,
    BaseAdapter,
)
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Import the module under test (will fail until implemented)
# =============================================================================

from workflow.adapters.busco import (
    BuscoAdapter,
    parse_short_summary,
    BUSCO_SUMMARY_PATTERN,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def adapter() -> BuscoAdapter:
    """Create a BuscoAdapter instance."""
    return BuscoAdapter()


@pytest.fixture
def basic_context(tmp_path: Path) -> AdapterContext:
    """Create a basic AdapterContext for testing."""
    # Create test input file
    proteins_file = tmp_path / "proteins.fa"
    proteins_file.write_text(">seq1\nMKTAYI\n>seq2\nMKTAYA\n")

    return AdapterContext(
        inputs={"proteins": proteins_file},
        outputs={
            "summary": tmp_path / "busco" / "short_summary.txt",
            "output_dir": tmp_path / "busco",
        },
        config={
            "busco": {
                "lineage": "primates_odb10",
                "mode": "proteins",
                "timeout": 1800,
            }
        },
        wildcards={"species": "mmur"},
        threads=8,
    )


@pytest.fixture
def sample_short_summary() -> str:
    """Sample BUSCO short_summary.txt content."""
    return """# BUSCO version is: 5.4.7
# The lineage dataset is: primates_odb10 (Creation date: 2024-01-08, number of BUSCOs: 13780)
# Summarized benchmarking in BUSCO notation for file proteins.fa
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:13780
        13118   Complete BUSCOs (C)
        12828   Complete and single-copy BUSCOs (S)
        290     Complete and duplicated BUSCOs (D)
        317     Fragmented BUSCOs (F)
        345     Missing BUSCOs (M)
        13780   Total BUSCO groups searched
"""


@pytest.fixture
def busco6_short_summary() -> str:
    """Sample BUSCO 6.x short_summary.txt content."""
    return """# BUSCO version is: 6.0.0
# The lineage dataset is: eukaryota_odb10 (Creation date: 2024-06-01, number of BUSCOs: 255)
# Summarized benchmarking in BUSCO notation for file proteins.fa
# BUSCO was run in mode: proteins

        ***** Results: *****

        C:98.4%[S:96.1%,D:2.3%],F:1.2%,M:0.4%,n:255
        251     Complete BUSCOs (C)
        245     Complete and single-copy BUSCOs (S)
        6       Complete and duplicated BUSCOs (D)
        3       Fragmented BUSCOs (F)
        1       Missing BUSCOs (M)
        255     Total BUSCO groups searched
"""


# =============================================================================
# Test BuscoAdapter.spec
# =============================================================================

class TestBuscoAdapterSpec:
    """Tests for BuscoAdapter.spec property."""

    def test_spec_returns_tool_spec(self, adapter: BuscoAdapter) -> None:
        """[P1] Given adapter, when accessing spec, then ToolSpec returned."""
        # When
        spec = adapter.spec

        # Then
        assert isinstance(spec, ToolSpec)

    def test_spec_name_is_busco(self, adapter: BuscoAdapter) -> None:
        """[P1] Given adapter, when accessing spec, then name is 'busco'."""
        # When
        spec = adapter.spec

        # Then
        assert spec.name == "busco"

    def test_spec_min_version(self, adapter: BuscoAdapter) -> None:
        """[P1] Given adapter, when accessing spec, then min_version is 5.0.0."""
        # When
        spec = adapter.spec

        # Then
        assert spec.min_version == "5.0.0"

    def test_spec_max_version(self, adapter: BuscoAdapter) -> None:
        """[P1] Given adapter, when accessing spec, then max_version is 6.99.99."""
        # When
        spec = adapter.spec

        # Then
        assert spec.max_version == "6.99.99"

    def test_spec_conda_env(self, adapter: BuscoAdapter) -> None:
        """[P1] Given adapter, when accessing spec, then conda_env is busco.yaml."""
        # When
        spec = adapter.spec

        # Then
        assert spec.conda_env == "busco.yaml"

    def test_adapter_is_base_adapter_subclass(self, adapter: BuscoAdapter) -> None:
        """[P1] Given adapter, then it is a BaseAdapter subclass."""
        # Then
        assert isinstance(adapter, BaseAdapter)


# =============================================================================
# Test BuscoAdapter.check_version
# =============================================================================

class TestBuscoAdapterCheckVersion:
    """Tests for BuscoAdapter.check_version method."""

    def test_check_version_5x(self, adapter: BuscoAdapter) -> None:
        """[P1] Given BUSCO 5.x installed, when check_version, then version returned."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                stdout="BUSCO 5.4.7\n",
                stderr="",
                returncode=0,
            )

            # When
            version = adapter.check_version()

            # Then
            assert version == "5.4.7"

    def test_check_version_6x(self, adapter: BuscoAdapter) -> None:
        """[P1] Given BUSCO 6.x installed, when check_version, then version returned."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                stdout="BUSCO 6.0.0\n",
                stderr="",
                returncode=0,
            )

            # When
            version = adapter.check_version()

            # Then
            assert version == "6.0.0"

    def test_check_version_from_stderr(self, adapter: BuscoAdapter) -> None:
        """[P2] Given version in stderr, when check_version, then version extracted."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                stdout="",
                stderr="BUSCO 5.5.0\n",
                returncode=0,
            )

            # When
            version = adapter.check_version()

            # Then
            assert version == "5.5.0"

    def test_check_version_unsupported_raises(self, adapter: BuscoAdapter) -> None:
        """[P1] Given unsupported version 4.x, when check_version, then E_TOOL_VERSION."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                stdout="BUSCO 4.1.4\n",
                stderr="",
                returncode=0,
            )

            # When/Then
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION

    def test_check_version_unsupported_7x_raises(self, adapter: BuscoAdapter) -> None:
        """[P2] Given unsupported version 7.x, when check_version, then E_TOOL_VERSION."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                stdout="BUSCO 7.0.0\n",
                stderr="",
                returncode=0,
            )

            # When/Then
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_VERSION

    def test_check_version_not_found_raises(self, adapter: BuscoAdapter) -> None:
        """[P1] Given BUSCO not installed, when check_version, then E_TOOL_NOT_FOUND."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = FileNotFoundError("busco not found")

            # When/Then
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND

    def test_check_version_parse_failure_raises(self, adapter: BuscoAdapter) -> None:
        """[P2] Given unparseable output, when check_version, then E_TOOL_NOT_FOUND."""
        # Given
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                stdout="Unknown output format",
                stderr="",
                returncode=0,
            )

            # When/Then
            with pytest.raises(CompGeneError) as exc_info:
                adapter.check_version()
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND


# =============================================================================
# Test BuscoAdapter.validate_inputs
# =============================================================================

class TestBuscoAdapterValidateInputs:
    """Tests for BuscoAdapter.validate_inputs method."""

    def test_validate_inputs_success(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given valid inputs, when validate_inputs, then no exception."""
        # When/Then (should not raise)
        adapter.validate_inputs(basic_context)

    def test_validate_inputs_missing_proteins_file(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P1] Given missing proteins file, when validate_inputs, then E_INPUT_MISSING."""
        # Given
        ctx = AdapterContext(
            inputs={"proteins": tmp_path / "nonexistent.fa"},
            outputs={},
            config={"busco": {"lineage": "primates_odb10"}},
            wildcards={},
            threads=1,
        )

        # When/Then
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING
        assert "protein" in exc_info.value.message.lower() or "not found" in exc_info.value.message.lower()

    def test_validate_inputs_missing_lineage_config(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P1] Given missing lineage config, when validate_inputs, then E_INPUT_MISSING."""
        # Given
        proteins_file = tmp_path / "proteins.fa"
        proteins_file.write_text(">seq\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins": proteins_file},
            outputs={},
            config={},  # No busco config
            wildcards={},
            threads=1,
        )

        # When/Then
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING
        assert "lineage" in exc_info.value.message.lower()

    def test_validate_inputs_empty_lineage(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P2] Given empty lineage string, when validate_inputs, then E_INPUT_MISSING."""
        # Given
        proteins_file = tmp_path / "proteins.fa"
        proteins_file.write_text(">seq\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins": proteins_file},
            outputs={},
            config={"busco": {"lineage": ""}},  # Empty lineage
            wildcards={},
            threads=1,
        )

        # When/Then
        with pytest.raises(CompGeneError) as exc_info:
            adapter.validate_inputs(ctx)
        assert exc_info.value.error_code == ErrorCode.E_INPUT_MISSING


# =============================================================================
# Test BuscoAdapter.build_command
# =============================================================================

class TestBuscoAdapterBuildCommand:
    """Tests for BuscoAdapter.build_command method."""

    def test_build_command_basic(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given basic config, when build_command, then correct command built."""
        # When
        cmd = adapter.build_command(basic_context)

        # Then
        assert cmd[0] == "busco"
        assert "-i" in cmd
        assert "-l" in cmd
        assert "primates_odb10" in cmd
        assert "-m" in cmd
        assert "proteins" in cmd
        assert "-c" in cmd
        assert "8" in cmd

    def test_build_command_includes_force_flag(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given config, when build_command, then -f flag included."""
        # When
        cmd = adapter.build_command(basic_context)

        # Then
        assert "-f" in cmd

    def test_build_command_with_download_path(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P2] Given download_path config, when build_command, then --download_path included."""
        # Given
        proteins_file = tmp_path / "proteins.fa"
        proteins_file.write_text(">seq\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins": proteins_file},
            outputs={"output_dir": tmp_path / "busco"},
            config={
                "busco": {
                    "lineage": "primates_odb10",
                    "download_path": "/data/busco_downloads",
                }
            },
            wildcards={"species": "mmur"},
            threads=4,
        )

        # When
        cmd = adapter.build_command(ctx)

        # Then
        assert "--download_path" in cmd
        assert "/data/busco_downloads" in cmd

    def test_build_command_with_offline_mode(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P2] Given offline=True, when build_command, then --offline included."""
        # Given
        proteins_file = tmp_path / "proteins.fa"
        proteins_file.write_text(">seq\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins": proteins_file},
            outputs={"output_dir": tmp_path / "busco"},
            config={
                "busco": {
                    "lineage": "primates_odb10",
                    "offline": True,
                }
            },
            wildcards={"species": "mmur"},
            threads=4,
        )

        # When
        cmd = adapter.build_command(ctx)

        # Then
        assert "--offline" in cmd

    def test_build_command_default_mode_is_proteins(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P2] Given no mode specified, when build_command, then mode is proteins."""
        # Given
        proteins_file = tmp_path / "proteins.fa"
        proteins_file.write_text(">seq\nMKT\n")

        ctx = AdapterContext(
            inputs={"proteins": proteins_file},
            outputs={"output_dir": tmp_path / "busco"},
            config={"busco": {"lineage": "primates_odb10"}},  # No mode
            wildcards={"species": "mmur"},
            threads=4,
        )

        # When
        cmd = adapter.build_command(ctx)

        # Then
        mode_idx = cmd.index("-m")
        assert cmd[mode_idx + 1] == "proteins"

    def test_build_command_uses_species_wildcard(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given species wildcard, when build_command, then -o uses species."""
        # When
        cmd = adapter.build_command(basic_context)

        # Then
        o_idx = cmd.index("-o")
        assert cmd[o_idx + 1] == "mmur"


# =============================================================================
# Test BuscoAdapter.expected_outputs
# =============================================================================

class TestBuscoAdapterExpectedOutputs:
    """Tests for BuscoAdapter.expected_outputs method."""

    def test_expected_outputs_returns_paths(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given context, when expected_outputs, then paths returned."""
        # When
        outputs = adapter.expected_outputs(basic_context)

        # Then
        assert isinstance(outputs, list)
        assert all(isinstance(p, Path) for p in outputs)

    def test_expected_outputs_includes_short_summary(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given context, when expected_outputs, then short_summary.txt included."""
        # When
        outputs = adapter.expected_outputs(basic_context)
        output_names = [p.name for p in outputs]

        # Then
        assert "short_summary.txt" in output_names


# =============================================================================
# Test BuscoAdapter.parse_outputs
# =============================================================================

class TestBuscoAdapterParseOutputs:
    """Tests for BuscoAdapter.parse_outputs method."""

    def test_parse_outputs_returns_run_result(
        self,
        adapter: BuscoAdapter,
        tmp_path: Path,
        sample_short_summary: str,
    ) -> None:
        """[P1] Given valid output, when parse_outputs, then RunResult returned."""
        # Given
        busco_dir = tmp_path / "busco"
        busco_dir.mkdir()
        summary_file = busco_dir / "short_summary.txt"
        summary_file.write_text(sample_short_summary)

        ctx = AdapterContext(
            inputs={},
            outputs={"summary": summary_file, "output_dir": busco_dir},
            config={},
            wildcards={"species": "mmur"},
            threads=1,
        )

        # When
        result = adapter.parse_outputs(ctx)

        # Then
        assert isinstance(result, RunResult)

    def test_parse_outputs_extracts_statistics(
        self,
        adapter: BuscoAdapter,
        tmp_path: Path,
        sample_short_summary: str,
    ) -> None:
        """[P1] Given valid output, when parse_outputs, then statistics in summary."""
        # Given
        busco_dir = tmp_path / "busco"
        busco_dir.mkdir()
        summary_file = busco_dir / "short_summary.txt"
        summary_file.write_text(sample_short_summary)

        ctx = AdapterContext(
            inputs={},
            outputs={"summary": summary_file, "output_dir": busco_dir},
            config={},
            wildcards={"species": "mmur"},
            threads=1,
        )

        # When
        result = adapter.parse_outputs(ctx)

        # Then
        assert "complete_pct" in result.summary
        assert result.summary["complete_pct"] == 95.2
        assert result.summary["single_copy_pct"] == 93.1
        assert result.summary["duplicated_pct"] == 2.1
        assert result.summary["fragmented_pct"] == 2.3
        assert result.summary["missing_pct"] == 2.5
        assert result.summary["total"] == 13780


# =============================================================================
# Test BuscoAdapter.timeout_seconds
# =============================================================================

class TestBuscoAdapterTimeoutSeconds:
    """Tests for BuscoAdapter.timeout_seconds method."""

    def test_timeout_seconds_default(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P1] Given no timeout config, when timeout_seconds, then default 1800."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},  # No busco config
            wildcards={},
            threads=1,
        )

        # When
        timeout = adapter.timeout_seconds(ctx)

        # Then
        assert timeout == 1800  # Default 30 minutes

    def test_timeout_seconds_from_config(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given timeout in config, when timeout_seconds, then config value used."""
        # Given - basic_context has timeout: 1800 in config

        # When
        timeout = adapter.timeout_seconds(basic_context)

        # Then
        assert timeout == 1800

    def test_timeout_seconds_custom(
        self, adapter: BuscoAdapter, tmp_path: Path
    ) -> None:
        """[P2] Given custom timeout, when timeout_seconds, then custom value used."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={"busco": {"timeout": 3600}},  # 1 hour
            wildcards={},
            threads=1,
        )

        # When
        timeout = adapter.timeout_seconds(ctx)

        # Then
        assert timeout == 3600


# =============================================================================
# Test BuscoAdapter.classify_error
# =============================================================================

class TestBuscoAdapterClassifyError:
    """Tests for BuscoAdapter.classify_error method."""

    def test_classify_error_lineage_not_found(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given lineage not found error, then E_INPUT_MISSING, not retryable."""
        # When
        code, retryable = adapter.classify_error(
            basic_context,
            1,
            "ERROR: The lineage primates_odb10 was not found in the downloads folder",
        )

        # Then
        assert code == ErrorCode.E_INPUT_MISSING
        assert retryable is False

    def test_classify_error_memory(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given memory error, then E_OOM, not retryable."""
        # When
        code, retryable = adapter.classify_error(
            basic_context,
            1,
            "MemoryError: Cannot allocate memory",
        )

        # Then
        assert code == ErrorCode.E_OOM
        assert retryable is False

    def test_classify_error_timeout_sigkill(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given SIGKILL (timeout), then E_TIMEOUT, retryable."""
        # When
        code, retryable = adapter.classify_error(
            basic_context,
            -9,  # SIGKILL
            "",
        )

        # Then
        assert code == ErrorCode.E_TIMEOUT
        assert retryable is True

    def test_classify_error_unknown(
        self, adapter: BuscoAdapter, basic_context: AdapterContext
    ) -> None:
        """[P1] Given unknown error, then E_NONZERO_EXIT, not retryable."""
        # When
        code, retryable = adapter.classify_error(
            basic_context,
            1,
            "Some unknown error",
        )

        # Then
        assert code == ErrorCode.E_NONZERO_EXIT
        assert retryable is False


# =============================================================================
# Test parse_short_summary helper function
# =============================================================================

class TestParseShortSummary:
    """Tests for parse_short_summary helper function."""

    def test_parse_short_summary_busco5(self, sample_short_summary: str) -> None:
        """[P1] Given BUSCO 5.x summary, when parsed, then correct values extracted."""
        # When
        result = parse_short_summary(sample_short_summary)

        # Then
        assert result["complete_pct"] == 95.2
        assert result["single_copy_pct"] == 93.1
        assert result["duplicated_pct"] == 2.1
        assert result["fragmented_pct"] == 2.3
        assert result["missing_pct"] == 2.5
        assert result["total"] == 13780

    def test_parse_short_summary_busco6(self, busco6_short_summary: str) -> None:
        """[P1] Given BUSCO 6.x summary, when parsed, then correct values extracted."""
        # When
        result = parse_short_summary(busco6_short_summary)

        # Then
        assert result["complete_pct"] == 98.4
        assert result["single_copy_pct"] == 96.1
        assert result["duplicated_pct"] == 2.3
        assert result["fragmented_pct"] == 1.2
        assert result["missing_pct"] == 0.4
        assert result["total"] == 255

    def test_parse_short_summary_invalid_raises(self) -> None:
        """[P1] Given invalid content, when parsed, then ValueError raised."""
        # When/Then
        with pytest.raises(ValueError) as exc_info:
            parse_short_summary("Invalid content without BUSCO results")
        assert "parse" in str(exc_info.value).lower() or "busco" in str(exc_info.value).lower()

    def test_busco_summary_pattern_matches(self) -> None:
        """[P2] Given BUSCO notation, when pattern applied, then matches."""
        # Given
        notation = "C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:13780"

        # When
        match = BUSCO_SUMMARY_PATTERN.search(notation)

        # Then
        assert match is not None
        assert match.group(1) == "95.2"
        assert match.group(2) == "93.1"
        assert match.group(3) == "2.1"
        assert match.group(4) == "2.3"
        assert match.group(5) == "2.5"
        assert match.group(6) == "13780"

    def test_parse_extracts_lineage(self, sample_short_summary: str) -> None:
        """[P2] Given summary with lineage, when parsed, then lineage extracted."""
        # When
        result = parse_short_summary(sample_short_summary)

        # Then
        assert result["lineage"] == "primates_odb10"

    def test_parse_extracts_version(self, sample_short_summary: str) -> None:
        """[P2] Given summary with version, when parsed, then version extracted."""
        # When
        result = parse_short_summary(sample_short_summary)

        # Then
        assert result["busco_version"] == "5.4.7"
