"""
Unit tests for workflow/adapters/base.py

Tests for ToolSpec, RunResult, AdapterContext, BaseAdapter, and error classification.
"""

import json
from pathlib import Path
from typing import Any

import pytest

from workflow.adapters.base import (
    ToolSpec,
    RunResult,
    AdapterContext,
    BaseAdapter,
    classify_common_errors,
)
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Test ToolSpec
# =============================================================================

class TestToolSpec:
    """Tests for ToolSpec dataclass."""

    def test_creates_with_required_fields(self) -> None:
        """[P1] Given required fields, when created, then all fields set."""
        # Given/When
        spec = ToolSpec(
            name="orthofinder",
            min_version="2.5.0"
        )

        # Then
        assert spec.name == "orthofinder"
        assert spec.min_version == "2.5.0"
        assert spec.max_version is None
        assert spec.conda_env == ""
        assert spec.description == ""

    def test_creates_with_all_fields(self) -> None:
        """[P1] Given all fields, when created, then all fields set."""
        # Given/When
        spec = ToolSpec(
            name="busco",
            min_version="5.0.0",
            max_version="6.99.99",
            conda_env="busco.yaml",
            description="BUSCO quality assessment"
        )

        # Then
        assert spec.name == "busco"
        assert spec.min_version == "5.0.0"
        assert spec.max_version == "6.99.99"
        assert spec.conda_env == "busco.yaml"
        assert spec.description == "BUSCO quality assessment"

    def test_version_compatible_within_range(self) -> None:
        """[P1] Given version in range, when check_version_compatible, then True."""
        # Given
        spec = ToolSpec(name="tool", min_version="2.0.0", max_version="2.99.99")

        # When/Then
        assert spec.check_version_compatible("2.0.0") is True
        assert spec.check_version_compatible("2.5.5") is True
        assert spec.check_version_compatible("2.99.99") is True

    def test_version_compatible_below_min(self) -> None:
        """[P1] Given version below min, when check_version_compatible, then False."""
        # Given
        spec = ToolSpec(name="tool", min_version="2.0.0")

        # When/Then
        assert spec.check_version_compatible("1.9.9") is False
        assert spec.check_version_compatible("1.0.0") is False

    def test_version_compatible_above_max(self) -> None:
        """[P1] Given version above max, when check_version_compatible, then False."""
        # Given
        spec = ToolSpec(name="tool", min_version="2.0.0", max_version="2.5.99")

        # When/Then
        assert spec.check_version_compatible("3.0.0") is False
        assert spec.check_version_compatible("2.6.0") is False

    def test_version_compatible_no_max(self) -> None:
        """[P2] Given no max version, when version above min, then True."""
        # Given
        spec = ToolSpec(name="tool", min_version="1.0.0")

        # When/Then
        assert spec.check_version_compatible("99.99.99") is True

    def test_parse_version_handles_various_formats(self) -> None:
        """[P2] Given various version formats, when _parse_version, then parsed correctly."""
        # Given
        spec = ToolSpec(name="tool", min_version="1.0.0")

        # When/Then
        assert spec._parse_version("2.5.5") == (2, 5, 5)
        assert spec._parse_version("1.0") == (1, 0)
        assert spec._parse_version("v2.1.3") == (2, 1, 3)
        assert spec._parse_version("version 1.2.3") == (1, 2, 3)

    def test_to_dict_includes_required_fields(self) -> None:
        """[P1] Given spec, when to_dict, then required fields included."""
        # Given
        spec = ToolSpec(name="tool", min_version="1.0.0", conda_env="tool.yaml")

        # When
        result = spec.to_dict()

        # Then
        assert result["name"] == "tool"
        assert result["min_version"] == "1.0.0"
        assert result["conda_env"] == "tool.yaml"

    def test_to_dict_excludes_empty_conda_env(self) -> None:
        """[P2] Given spec without conda_env, when to_dict, then conda_env excluded."""
        # Given
        spec = ToolSpec(name="tool", min_version="1.0.0")

        # When
        result = spec.to_dict()

        # Then
        assert "conda_env" not in result

    def test_to_dict_excludes_none_max_version(self) -> None:
        """[P2] Given spec without max_version, when to_dict, then max_version excluded."""
        # Given
        spec = ToolSpec(name="tool", min_version="1.0.0")

        # When
        result = spec.to_dict()

        # Then
        assert "max_version" not in result

    def test_to_dict_includes_max_version_when_set(self) -> None:
        """[P2] Given spec with max_version, when to_dict, then max_version included."""
        # Given
        spec = ToolSpec(name="tool", min_version="1.0.0", max_version="2.0.0")

        # When
        result = spec.to_dict()

        # Then
        assert result["max_version"] == "2.0.0"


# =============================================================================
# Test RunResult
# =============================================================================

class TestRunResult:
    """Tests for RunResult dataclass."""

    def test_creates_with_required_fields(self) -> None:
        """[P1] Given required fields, when created, then all fields set."""
        # Given/When
        result = RunResult(
            outputs={"output1": Path("results/output.tsv")},
            summary={"count": 100}
        )

        # Then
        assert "output1" in result.outputs
        assert result.outputs["output1"] == Path("results/output.tsv")
        assert result.summary["count"] == 100
        assert result.metrics is None
        assert result.warnings is None

    def test_creates_with_all_fields(self) -> None:
        """[P1] Given all fields, when created, then all fields set."""
        # Given/When
        result = RunResult(
            outputs={"orthogroups": Path("results/og.tsv")},
            summary={"total": 1000, "species": 5},
            metrics={"runtime_seconds": 3600.5, "memory_mb": 8000},
            warnings=["Low memory warning"]
        )

        # Then
        assert result.outputs["orthogroups"] == Path("results/og.tsv")
        assert result.summary["total"] == 1000
        assert result.metrics["runtime_seconds"] == 3600.5
        assert result.warnings == ["Low memory warning"]

    def test_to_dict_converts_paths_to_strings(self) -> None:
        """[P1] Given result with paths, when to_dict, then paths are strings."""
        # Given
        result = RunResult(
            outputs={"out": Path("/path/to/output.tsv")},
            summary={"count": 1}
        )

        # When
        d = result.to_dict()

        # Then
        assert d["outputs"]["out"] == "/path/to/output.tsv"
        assert isinstance(d["outputs"]["out"], str)

    def test_to_dict_excludes_none_optionals(self) -> None:
        """[P2] Given result without optionals, when to_dict, then None excluded."""
        # Given
        result = RunResult(outputs={}, summary={})

        # When
        d = result.to_dict()

        # Then
        assert "metrics" not in d
        assert "warnings" not in d

    def test_to_dict_includes_optionals_when_set(self) -> None:
        """[P2] Given result with optionals, when to_dict, then optionals included."""
        # Given
        result = RunResult(
            outputs={},
            summary={},
            metrics={"time": 10},
            warnings=["warn1"]
        )

        # When
        d = result.to_dict()

        # Then
        assert d["metrics"] == {"time": 10}
        assert d["warnings"] == ["warn1"]

    def test_to_json_returns_valid_json(self) -> None:
        """[P1] Given result, when to_json, then valid JSON returned."""
        # Given
        result = RunResult(
            outputs={"file": Path("output.txt")},
            summary={"count": 42}
        )

        # When
        json_str = result.to_json()

        # Then
        parsed = json.loads(json_str)
        assert parsed["outputs"]["file"] == "output.txt"
        assert parsed["summary"]["count"] == 42


# =============================================================================
# Test AdapterContext
# =============================================================================

class TestAdapterContext:
    """Tests for AdapterContext dataclass."""

    def test_creates_with_required_fields(self) -> None:
        """[P1] Given required fields, when created, then all fields set."""
        # Given/When
        ctx = AdapterContext(
            inputs={"input1": Path("data/input.fa")},
            outputs={"output1": Path("results/output.tsv")},
            config={"threads": 8},
            wildcards={"species": "mmur"},
            threads=8
        )

        # Then
        assert ctx.inputs["input1"] == Path("data/input.fa")
        assert ctx.outputs["output1"] == Path("results/output.tsv")
        assert ctx.config["threads"] == 8
        assert ctx.wildcards["species"] == "mmur"
        assert ctx.threads == 8
        assert ctx.logger is None
        assert ctx.temp_dir is None

    def test_get_input_returns_path(self) -> None:
        """[P1] Given valid input name, when get_input, then path returned."""
        # Given
        ctx = AdapterContext(
            inputs={"proteins": Path("data/proteins.fa")},
            outputs={},
            config={},
            wildcards={},
            threads=1
        )

        # When
        path = ctx.get_input("proteins")

        # Then
        assert path == Path("data/proteins.fa")

    def test_get_input_raises_for_missing(self) -> None:
        """[P1] Given invalid input name, when get_input, then KeyError raised."""
        # Given
        ctx = AdapterContext(
            inputs={"existing": Path("data/existing.fa")},
            outputs={},
            config={},
            wildcards={},
            threads=1
        )

        # When/Then
        with pytest.raises(KeyError) as exc_info:
            ctx.get_input("nonexistent")
        assert "nonexistent" in str(exc_info.value)

    def test_get_output_returns_path(self) -> None:
        """[P1] Given valid output name, when get_output, then path returned."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={"results": Path("results/output.tsv")},
            config={},
            wildcards={},
            threads=1
        )

        # When
        path = ctx.get_output("results")

        # Then
        assert path == Path("results/output.tsv")

    def test_get_output_raises_for_missing(self) -> None:
        """[P1] Given invalid output name, when get_output, then KeyError raised."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={"existing": Path("results/existing.tsv")},
            config={},
            wildcards={},
            threads=1
        )

        # When/Then
        with pytest.raises(KeyError) as exc_info:
            ctx.get_output("nonexistent")
        assert "nonexistent" in str(exc_info.value)

    def test_get_config_returns_value(self) -> None:
        """[P2] Given config key, when get_config, then value returned."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={"memory_mb": 16000},
            wildcards={},
            threads=1
        )

        # When
        value = ctx.get_config("memory_mb")

        # Then
        assert value == 16000

    def test_get_config_returns_default_for_missing(self) -> None:
        """[P2] Given missing key, when get_config with default, then default returned."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},
            wildcards={},
            threads=1
        )

        # When
        value = ctx.get_config("missing_key", default=100)

        # Then
        assert value == 100

    def test_get_wildcard_returns_value(self) -> None:
        """[P2] Given wildcard name, when get_wildcard, then value returned."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},
            wildcards={"species_set": "lemur_macaque"},
            threads=1
        )

        # When
        value = ctx.get_wildcard("species_set")

        # Then
        assert value == "lemur_macaque"

    def test_get_wildcard_raises_for_missing(self) -> None:
        """[P2] Given missing wildcard, when get_wildcard, then KeyError raised."""
        # Given
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},
            wildcards={"existing": "value"},
            threads=1
        )

        # When/Then
        with pytest.raises(KeyError) as exc_info:
            ctx.get_wildcard("nonexistent")
        assert "nonexistent" in str(exc_info.value)


# =============================================================================
# Test BaseAdapter (via MockAdapter)
# =============================================================================

class MockAdapter(BaseAdapter):
    """Mock adapter for testing BaseAdapter interface."""

    def __init__(self, tool_spec: ToolSpec | None = None):
        self._spec = tool_spec or ToolSpec(
            name="mock_tool",
            min_version="1.0.0",
            max_version="2.0.0",
            conda_env="mock.yaml",
            description="Mock tool for testing"
        )
        self._version = "1.5.0"

    @property
    def spec(self) -> ToolSpec:
        return self._spec

    def check_version(self) -> str:
        return self._version

    def validate_inputs(self, ctx: AdapterContext) -> None:
        # Check all inputs exist
        for name, path in ctx.inputs.items():
            if not path.exists():
                raise CompGeneError(
                    ErrorCode.E_INPUT_MISSING,
                    f"Input '{name}' not found",
                    details=str(path)
                )

    def build_command(self, ctx: AdapterContext) -> list[str]:
        cmd = ["mock_tool"]
        for name, path in ctx.inputs.items():
            cmd.extend(["-i", str(path)])
        for name, path in ctx.outputs.items():
            cmd.extend(["-o", str(path)])
        cmd.extend(["-t", str(ctx.threads)])
        return cmd

    def expected_outputs(self, ctx: AdapterContext) -> list[Path]:
        return list(ctx.outputs.values())

    def parse_outputs(self, ctx: AdapterContext) -> RunResult:
        return RunResult(
            outputs=ctx.outputs,
            summary={"status": "success"},
            metrics={"threads": ctx.threads}
        )

    def timeout_seconds(self, ctx: AdapterContext) -> int:
        # Base timeout + time per input
        return 300 + (len(ctx.inputs) * 60)

    def classify_error(
        self,
        ctx: AdapterContext,
        returncode: int,
        stderr: str
    ) -> tuple[ErrorCode, bool]:
        return classify_common_errors(returncode, stderr)


class TestBaseAdapter:
    """Tests for BaseAdapter abstract class via MockAdapter."""

    def test_can_instantiate_subclass(self) -> None:
        """[P1] Given subclass implementing all methods, when instantiated, then works."""
        # Given/When
        adapter = MockAdapter()

        # Then
        assert adapter is not None
        assert isinstance(adapter, BaseAdapter)

    def test_spec_returns_tool_spec(self) -> None:
        """[P1] Given adapter, when accessing spec, then ToolSpec returned."""
        # Given
        adapter = MockAdapter()

        # When
        spec = adapter.spec

        # Then
        assert isinstance(spec, ToolSpec)
        assert spec.name == "mock_tool"

    def test_check_version_returns_string(self) -> None:
        """[P1] Given adapter, when check_version, then version string returned."""
        # Given
        adapter = MockAdapter()

        # When
        version = adapter.check_version()

        # Then
        assert version == "1.5.0"

    def test_build_command_returns_list(self) -> None:
        """[P1] Given context, when build_command, then command list returned."""
        # Given
        adapter = MockAdapter()
        ctx = AdapterContext(
            inputs={"input": Path("data/input.fa")},
            outputs={"output": Path("results/output.tsv")},
            config={},
            wildcards={},
            threads=4
        )

        # When
        cmd = adapter.build_command(ctx)

        # Then
        assert isinstance(cmd, list)
        assert cmd[0] == "mock_tool"
        assert "-i" in cmd
        assert "-o" in cmd
        assert "-t" in cmd
        assert "4" in cmd

    def test_expected_outputs_returns_paths(self) -> None:
        """[P1] Given context, when expected_outputs, then path list returned."""
        # Given
        adapter = MockAdapter()
        ctx = AdapterContext(
            inputs={},
            outputs={"out1": Path("results/out1.tsv"), "out2": Path("results/out2.tsv")},
            config={},
            wildcards={},
            threads=1
        )

        # When
        outputs = adapter.expected_outputs(ctx)

        # Then
        assert len(outputs) == 2
        assert all(isinstance(p, Path) for p in outputs)

    def test_parse_outputs_returns_run_result(self) -> None:
        """[P1] Given context, when parse_outputs, then RunResult returned."""
        # Given
        adapter = MockAdapter()
        ctx = AdapterContext(
            inputs={},
            outputs={"out": Path("results/out.tsv")},
            config={},
            wildcards={},
            threads=8
        )

        # When
        result = adapter.parse_outputs(ctx)

        # Then
        assert isinstance(result, RunResult)
        assert result.summary["status"] == "success"
        assert result.metrics["threads"] == 8

    def test_timeout_seconds_returns_int(self) -> None:
        """[P1] Given context, when timeout_seconds, then int returned."""
        # Given
        adapter = MockAdapter()
        ctx = AdapterContext(
            inputs={"in1": Path("a"), "in2": Path("b")},
            outputs={},
            config={},
            wildcards={},
            threads=1
        )

        # When
        timeout = adapter.timeout_seconds(ctx)

        # Then
        assert isinstance(timeout, int)
        assert timeout == 300 + (2 * 60)  # base + 2 inputs

    def test_classify_error_returns_tuple(self) -> None:
        """[P1] Given error info, when classify_error, then (ErrorCode, bool) returned."""
        # Given
        adapter = MockAdapter()
        ctx = AdapterContext(
            inputs={},
            outputs={},
            config={},
            wildcards={},
            threads=1
        )

        # When
        code, retryable = adapter.classify_error(ctx, 1, "some error")

        # Then
        assert isinstance(code, ErrorCode)
        assert isinstance(retryable, bool)


# =============================================================================
# Test classify_common_errors
# =============================================================================

class TestClassifyCommonErrors:
    """Tests for classify_common_errors helper function."""

    def test_tool_not_found_returncode_127(self) -> None:
        """[P1] Given returncode 127, then E_TOOL_NOT_FOUND, not retryable."""
        code, retryable = classify_common_errors(127, "")
        assert code == ErrorCode.E_TOOL_NOT_FOUND
        assert retryable is False

    def test_permission_denied_returncode_126(self) -> None:
        """[P2] Given returncode 126, then E_TOOL_NOT_FOUND, not retryable."""
        code, retryable = classify_common_errors(126, "")
        assert code == ErrorCode.E_TOOL_NOT_FOUND
        assert retryable is False

    def test_timeout_in_stderr(self) -> None:
        """[P1] Given timeout in stderr, then E_TIMEOUT, retryable."""
        code, retryable = classify_common_errors(1, "Process timed out after 30 minutes")
        assert code == ErrorCode.E_TIMEOUT
        assert retryable is True

    def test_out_of_memory(self) -> None:
        """[P1] Given OOM message, then E_OOM, not retryable."""
        code, retryable = classify_common_errors(1, "Out of memory error")
        assert code == ErrorCode.E_OOM
        assert retryable is False

    def test_cannot_allocate_memory(self) -> None:
        """[P2] Given cannot allocate message, then E_OOM, not retryable."""
        code, retryable = classify_common_errors(1, "Cannot allocate memory")
        assert code == ErrorCode.E_OOM
        assert retryable is False

    def test_bad_alloc(self) -> None:
        """[P2] Given bad_alloc message, then E_OOM, not retryable."""
        code, retryable = classify_common_errors(1, "std::bad_alloc")
        assert code == ErrorCode.E_OOM
        assert retryable is False

    def test_disk_full(self) -> None:
        """[P1] Given disk full message, then E_DISK_FULL, not retryable."""
        code, retryable = classify_common_errors(1, "No space left on device")
        assert code == ErrorCode.E_DISK_FULL
        assert retryable is False

    def test_rate_limit(self) -> None:
        """[P1] Given rate limit message, then E_NET_RATE_LIMIT, retryable."""
        code, retryable = classify_common_errors(1, "API rate limit exceeded")
        assert code == ErrorCode.E_NET_RATE_LIMIT
        assert retryable is True

    def test_connection_refused(self) -> None:
        """[P2] Given connection refused, then E_NET_RATE_LIMIT, retryable."""
        code, retryable = classify_common_errors(1, "Connection refused")
        assert code == ErrorCode.E_NET_RATE_LIMIT
        assert retryable is True

    def test_input_format_error(self) -> None:
        """[P1] Given format error message, then E_INPUT_FORMAT, not retryable."""
        code, retryable = classify_common_errors(1, "Invalid file format error")
        assert code == ErrorCode.E_INPUT_FORMAT
        assert retryable is False

    def test_file_not_found(self) -> None:
        """[P1] Given file not found message, then E_INPUT_MISSING, not retryable."""
        code, retryable = classify_common_errors(1, "File not found: input.fa")
        assert code == ErrorCode.E_INPUT_MISSING
        assert retryable is False

    def test_no_such_file(self) -> None:
        """[P2] Given no such file message, then E_INPUT_MISSING, not retryable."""
        code, retryable = classify_common_errors(1, "No such file or directory")
        assert code == ErrorCode.E_INPUT_MISSING
        assert retryable is False

    def test_default_nonzero_exit(self) -> None:
        """[P1] Given unknown error, then E_NONZERO_EXIT, not retryable."""
        code, retryable = classify_common_errors(1, "Unknown error occurred")
        assert code == ErrorCode.E_NONZERO_EXIT
        assert retryable is False
