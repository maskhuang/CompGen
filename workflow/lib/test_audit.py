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
    get_checkpoint_path,
    mark_rule_complete,
    is_rule_complete,
    get_rule_completion_time,
    clear_checkpoint,
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


# =============================================================================
# Test Checkpoint Functions (Story 1.7)
# =============================================================================

class TestGetCheckpointPath:
    """Tests for get_checkpoint_path function."""

    def test_generates_correct_path(self, tmp_path: Path) -> None:
        """[P1] Given rule and wildcards, then correct .done path returned."""
        # Given
        rule = "orthofinder"
        wildcards = {"species_set": "lemur"}

        # When
        path = get_checkpoint_path(rule, wildcards, tmp_path)

        # Then
        assert path.suffix == ".done"
        assert "orthofinder" in str(path)
        assert "species_set=lemur" in str(path)

    def test_uses_default_for_empty_wildcards(self, tmp_path: Path) -> None:
        """[P2] Given empty wildcards, then 'default' in filename."""
        # Given
        rule = "test_rule"
        wildcards: dict[str, str] = {}

        # When
        path = get_checkpoint_path(rule, wildcards, tmp_path)

        # Then
        assert "default.done" in str(path)

    def test_sorts_wildcards_deterministically(self, tmp_path: Path) -> None:
        """[P2] Given multiple wildcards, then sorted alphabetically in path."""
        # Given
        rule = "test"
        wildcards = {"z_last": "1", "a_first": "2", "m_middle": "3"}

        # When
        path = get_checkpoint_path(rule, wildcards, tmp_path)

        # Then
        path_str = str(path)
        # Should be a_first before m_middle before z_last
        assert path_str.index("a_first") < path_str.index("m_middle")
        assert path_str.index("m_middle") < path_str.index("z_last")


class TestMarkRuleComplete:
    """Tests for mark_rule_complete function."""

    def test_creates_checkpoint_file(self, tmp_path: Path) -> None:
        """[P1] Given rule info, when marked complete, then .done file created."""
        # Given
        rule = "test_rule"
        wildcards = {"species": "mmur"}

        # When
        path = mark_rule_complete(rule, wildcards, tmp_path)

        # Then
        assert path.exists()
        assert path.suffix == ".done"

    def test_writes_timestamp(self, tmp_path: Path) -> None:
        """[P1] Given call, when marked complete, then timestamp in file."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}

        # When
        path = mark_rule_complete(rule, wildcards, tmp_path)

        # Then
        content = path.read_text().strip()
        assert "T" in content  # ISO format contains T
        assert "+" in content or "Z" in content  # Timezone info

    def test_uses_provided_timestamp(self, tmp_path: Path) -> None:
        """[P2] Given explicit timestamp, when marked, then that timestamp used."""
        # Given
        rule = "test"
        wildcards = {}
        timestamp = "2026-01-28T12:00:00+00:00"

        # When
        path = mark_rule_complete(rule, wildcards, tmp_path, timestamp=timestamp)

        # Then
        content = path.read_text().strip()
        assert content == timestamp

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        """[P2] Given nested path, when marked, then directories created."""
        # Given
        rule = "nested_rule"
        wildcards = {"key": "value"}

        # When
        path = mark_rule_complete(rule, wildcards, tmp_path / "deep" / "path")

        # Then
        assert path.exists()


class TestIsRuleComplete:
    """Tests for is_rule_complete function."""

    def test_returns_true_when_complete(self, tmp_path: Path) -> None:
        """[P1] Given checkpoint exists, when checked, then True returned."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        mark_rule_complete(rule, wildcards, tmp_path)

        # When
        result = is_rule_complete(rule, wildcards, tmp_path)

        # Then
        assert result is True

    def test_returns_false_when_not_complete(self, tmp_path: Path) -> None:
        """[P1] Given no checkpoint, when checked, then False returned."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        # Not marked complete

        # When
        result = is_rule_complete(rule, wildcards, tmp_path)

        # Then
        assert result is False

    def test_distinguishes_different_wildcards(self, tmp_path: Path) -> None:
        """[P2] Given same rule different wildcards, then separate completion status."""
        # Given
        rule = "test"
        mark_rule_complete(rule, {"species": "mmur"}, tmp_path)
        # Not marked for lcat

        # When/Then
        assert is_rule_complete(rule, {"species": "mmur"}, tmp_path) is True
        assert is_rule_complete(rule, {"species": "lcat"}, tmp_path) is False


class TestGetRuleCompletionTime:
    """Tests for get_rule_completion_time function."""

    def test_returns_timestamp_when_complete(self, tmp_path: Path) -> None:
        """[P1] Given checkpoint exists, when called, then timestamp returned."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        timestamp = "2026-01-28T12:00:00+00:00"
        mark_rule_complete(rule, wildcards, tmp_path, timestamp=timestamp)

        # When
        result = get_rule_completion_time(rule, wildcards, tmp_path)

        # Then
        assert result == timestamp

    def test_returns_none_when_not_complete(self, tmp_path: Path) -> None:
        """[P1] Given no checkpoint, when called, then None returned."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        # Not marked complete

        # When
        result = get_rule_completion_time(rule, wildcards, tmp_path)

        # Then
        assert result is None


class TestClearCheckpoint:
    """Tests for clear_checkpoint function."""

    def test_removes_checkpoint_file(self, tmp_path: Path) -> None:
        """[P1] Given checkpoint exists, when cleared, then file removed."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        path = mark_rule_complete(rule, wildcards, tmp_path)
        assert path.exists()

        # When
        result = clear_checkpoint(rule, wildcards, tmp_path)

        # Then
        assert result is True
        assert not path.exists()

    def test_returns_false_when_no_checkpoint(self, tmp_path: Path) -> None:
        """[P1] Given no checkpoint, when cleared, then False returned."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        # Not marked complete

        # When
        result = clear_checkpoint(rule, wildcards, tmp_path)

        # Then
        assert result is False

    def test_allows_re_completion_after_clear(self, tmp_path: Path) -> None:
        """[P2] Given cleared checkpoint, when marked again, then new checkpoint created."""
        # Given
        rule = "test"
        wildcards = {"key": "value"}
        mark_rule_complete(rule, wildcards, tmp_path, timestamp="2026-01-01T00:00:00Z")
        clear_checkpoint(rule, wildcards, tmp_path)

        # When
        new_path = mark_rule_complete(rule, wildcards, tmp_path, timestamp="2026-01-02T00:00:00Z")

        # Then
        assert new_path.exists()
        assert get_rule_completion_time(rule, wildcards, tmp_path) == "2026-01-02T00:00:00Z"
