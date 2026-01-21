"""
Unit tests for workflow/lib/audit.py

Tests audit metadata collection and writing.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from workflow.lib.audit import (
    RunMetadata,
    get_audit_path,
    collect_run_metadata,
    write_run_json,
    create_and_write_audit,
)


# =============================================================================
# Test RunMetadata
# =============================================================================

class TestRunMetadata:
    """Tests for RunMetadata dataclass."""

    def test_creates_with_required_fields(self) -> None:
        """[P1] Given required fields, when created, then all fields set."""
        # Given/When
        metadata = RunMetadata(
            rule="test_rule",
            wildcards={"species": "mmur"},
            cmd=["echo", "hello"],
            tool_version="1.0.0",
            input_checksums={"input.txt": "sha256:abc123"},
            threads=4,
            runtime_seconds=10.5,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00"
        )

        # Then
        assert metadata.rule == "test_rule"
        assert metadata.wildcards == {"species": "mmur"}
        assert metadata.threads == 4
        assert metadata.exit_code == 0

    def test_optional_fields_default_to_none(self) -> None:
        """[P1] Given only required fields, then optional fields are None."""
        # Given/When
        metadata = RunMetadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_checksums={},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00"
        )

        # Then
        assert metadata.output_checksums is None
        assert metadata.error_code is None
        assert metadata.error_message is None

    def test_to_dict_excludes_none_values(self) -> None:
        """[P1] Given metadata with None fields, when to_dict called, then Nones excluded."""
        # Given
        metadata = RunMetadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_checksums={},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00",
            output_checksums=None,  # Should be excluded
            error_code=None         # Should be excluded
        )

        # When
        result = metadata.to_dict()

        # Then
        assert "output_checksums" not in result
        assert "error_code" not in result
        assert "error_message" not in result

    def test_to_dict_includes_non_none_optionals(self) -> None:
        """[P2] Given metadata with optional fields set, when to_dict called, then included."""
        # Given
        metadata = RunMetadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_checksums={},
            threads=1,
            runtime_seconds=0.0,
            exit_code=1,
            timestamp="2026-01-20T10:00:00+00:00",
            error_code="E_TIMEOUT",
            error_message="Process timed out"
        )

        # When
        result = metadata.to_dict()

        # Then
        assert result["error_code"] == "E_TIMEOUT"
        assert result["error_message"] == "Process timed out"

    def test_to_json_returns_valid_json(self) -> None:
        """[P1] Given metadata, when to_json called, then valid JSON returned."""
        # Given
        metadata = RunMetadata(
            rule="test",
            wildcards={"species": "mmur"},
            cmd=["echo", "test"],
            tool_version="1.0",
            input_checksums={"file": "sha256:abc"},
            threads=4,
            runtime_seconds=10.5,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00"
        )

        # When
        json_str = metadata.to_json()

        # Then
        parsed = json.loads(json_str)
        assert parsed["rule"] == "test"
        assert parsed["threads"] == 4


# =============================================================================
# Test get_audit_path
# =============================================================================

class TestGetAuditPath:
    """Tests for get_audit_path function."""

    def test_generates_correct_path(self, tmp_path: Path) -> None:
        """[P1] Given rule and wildcards, then correct path returned."""
        # Given
        rule = "orthofinder"
        wildcards = {"species_set": "lemur"}

        # When
        path = get_audit_path(rule, wildcards, tmp_path)

        # Then
        assert path.suffix == ".json"
        assert "orthofinder" in str(path)
        assert "species_set=lemur" in str(path)

    def test_uses_default_for_empty_wildcards(self, tmp_path: Path) -> None:
        """[P2] Given empty wildcards, then 'default' in filename."""
        # Given
        rule = "test_rule"
        wildcards: dict[str, str] = {}

        # When
        path = get_audit_path(rule, wildcards, tmp_path)

        # Then
        assert "default" in path.stem

    def test_path_ends_with_run_json(self, tmp_path: Path) -> None:
        """[P2] Given any input, then path ends with .run.json."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}

        # When
        path = get_audit_path(rule, wildcards, tmp_path)

        # Then
        assert str(path).endswith(".run.json")


# =============================================================================
# Test collect_run_metadata
# =============================================================================

class TestCollectRunMetadata:
    """Tests for collect_run_metadata function."""

    def test_collects_basic_metadata(self, tmp_path: Path) -> None:
        """[P1] Given run info, when collected, then all fields populated."""
        # Given
        input_file = tmp_path / "input.txt"
        input_file.write_text("test content")

        # When
        metadata = collect_run_metadata(
            rule="test_rule",
            wildcards={"species": "mmur"},
            cmd=["tool", "-i", str(input_file)],
            tool_version="2.0.0",
            input_paths={"input": input_file},
            threads=8,
            runtime_seconds=120.5,
            exit_code=0
        )

        # Then
        assert metadata.rule == "test_rule"
        assert metadata.wildcards == {"species": "mmur"}
        assert metadata.cmd == ["tool", "-i", str(input_file)]
        assert metadata.tool_version == "2.0.0"
        assert metadata.threads == 8
        assert metadata.runtime_seconds == 120.5
        assert metadata.exit_code == 0

    def test_computes_input_checksums(self, tmp_path: Path) -> None:
        """[P1] Given input paths, when collected, then checksums computed."""
        # Given
        input_file = tmp_path / "input.txt"
        input_file.write_text("test content")

        # When
        metadata = collect_run_metadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_paths={"input": input_file},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0
        )

        # Then
        assert "input" in metadata.input_checksums
        assert metadata.input_checksums["input"].startswith("sha256:")

    def test_generates_timestamp(self, tmp_path: Path) -> None:
        """[P1] Given call, when collected, then timestamp generated."""
        # Given
        input_file = tmp_path / "input.txt"
        input_file.write_text("test")

        # When
        metadata = collect_run_metadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_paths={"input": input_file},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0
        )

        # Then
        assert metadata.timestamp is not None
        assert "T" in metadata.timestamp  # ISO format

    def test_includes_output_checksums_when_provided(self, tmp_path: Path) -> None:
        """[P2] Given output paths, when collected, then output checksums computed."""
        # Given
        input_file = tmp_path / "input.txt"
        output_file = tmp_path / "output.txt"
        input_file.write_text("input")
        output_file.write_text("output")

        # When
        metadata = collect_run_metadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_paths={"input": input_file},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0,
            output_paths={"output": output_file}
        )

        # Then
        assert metadata.output_checksums is not None
        assert "output" in metadata.output_checksums

    def test_includes_error_info_when_provided(self, tmp_path: Path) -> None:
        """[P2] Given error info, when collected, then error fields populated."""
        # Given
        input_file = tmp_path / "input.txt"
        input_file.write_text("test")

        # When
        metadata = collect_run_metadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_paths={"input": input_file},
            threads=1,
            runtime_seconds=0.0,
            exit_code=1,
            error_code="E_TIMEOUT",
            error_message="Process timed out after 30s"
        )

        # Then
        assert metadata.error_code == "E_TIMEOUT"
        assert metadata.error_message == "Process timed out after 30s"


# =============================================================================
# Test write_run_json
# =============================================================================

class TestWriteRunJson:
    """Tests for write_run_json function."""

    def test_writes_json_file(self, tmp_path: Path) -> None:
        """[P1] Given metadata, when written, then valid JSON file created."""
        # Given
        metadata = RunMetadata(
            rule="test_rule",
            wildcards={"species": "mmur"},
            cmd=["echo", "test"],
            tool_version="1.0",
            input_checksums={"input": "sha256:abc"},
            threads=4,
            runtime_seconds=10.0,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00"
        )

        # When
        path = write_run_json(metadata, tmp_path)

        # Then
        assert path.exists()
        content = json.loads(path.read_text())
        assert content["rule"] == "test_rule"

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        """[P1] Given nested path, when written, then directories created."""
        # Given
        metadata = RunMetadata(
            rule="nested_rule",
            wildcards={"key": "value"},
            cmd=["cmd"],
            tool_version="1.0",
            input_checksums={},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00"
        )

        # When
        path = write_run_json(metadata, tmp_path / "deep" / "nested")

        # Then
        assert path.exists()

    def test_returns_written_path(self, tmp_path: Path) -> None:
        """[P2] Given metadata, when written, then path returned."""
        # Given
        metadata = RunMetadata(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_checksums={},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0,
            timestamp="2026-01-20T10:00:00+00:00"
        )

        # When
        path = write_run_json(metadata, tmp_path)

        # Then
        assert path.suffix == ".json"
        assert "test" in str(path)


# =============================================================================
# Test create_and_write_audit
# =============================================================================

class TestCreateAndWriteAudit:
    """Tests for create_and_write_audit convenience function."""

    def test_creates_and_writes_in_one_step(self, tmp_path: Path) -> None:
        """[P1] Given run info, when called, then metadata created and written."""
        # Given
        input_file = tmp_path / "input.txt"
        input_file.write_text("test content")

        # When
        metadata, path = create_and_write_audit(
            rule="combined_test",
            wildcards={"species": "mmur"},
            cmd=["tool", "arg"],
            tool_version="2.0",
            input_paths={"input": input_file},
            threads=4,
            runtime_seconds=60.0,
            exit_code=0,
            meta_dir=tmp_path
        )

        # Then
        assert metadata.rule == "combined_test"
        assert path.exists()
        content = json.loads(path.read_text())
        assert content["rule"] == "combined_test"

    def test_returns_both_metadata_and_path(self, tmp_path: Path) -> None:
        """[P2] Given call, then tuple of (RunMetadata, Path) returned."""
        # Given
        input_file = tmp_path / "input.txt"
        input_file.write_text("test")

        # When
        result = create_and_write_audit(
            rule="test",
            wildcards={},
            cmd=["cmd"],
            tool_version="1.0",
            input_paths={"input": input_file},
            threads=1,
            runtime_seconds=0.0,
            exit_code=0,
            meta_dir=tmp_path
        )

        # Then
        assert len(result) == 2
        assert isinstance(result[0], RunMetadata)
        assert isinstance(result[1], Path)
