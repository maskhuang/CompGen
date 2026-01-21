"""
Unit tests for workflow/lib/io.py

Tests atomic file writing functionality.
"""

import json
import tempfile
from pathlib import Path

import pytest

from workflow.lib.io import atomic_write, atomic_write_json, atomic_append


# =============================================================================
# Test atomic_write
# =============================================================================

class TestAtomicWrite:
    """Tests for atomic_write function."""

    def test_writes_content_to_file(self, tmp_path: Path) -> None:
        """[P1] Given content, when atomic_write called, then file contains content."""
        # Given
        file_path = tmp_path / "test.txt"
        content = "Hello, World!"

        # When
        atomic_write(file_path, content)

        # Then
        assert file_path.exists()
        assert file_path.read_text() == content

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        """[P1] Given nested path, when atomic_write called, then directories created."""
        # Given
        file_path = tmp_path / "subdir" / "nested" / "test.txt"
        content = "Nested content"

        # When
        atomic_write(file_path, content)

        # Then
        assert file_path.exists()
        assert file_path.read_text() == content

    def test_overwrites_existing_file(self, tmp_path: Path) -> None:
        """[P1] Given existing file, when atomic_write called, then file overwritten."""
        # Given
        file_path = tmp_path / "test.txt"
        file_path.write_text("Old content")
        new_content = "New content"

        # When
        atomic_write(file_path, new_content)

        # Then
        assert file_path.read_text() == new_content

    def test_no_temp_file_remains_after_success(self, tmp_path: Path) -> None:
        """[P2] Given successful write, when complete, then no .tmp file remains."""
        # Given
        file_path = tmp_path / "test.txt"
        temp_path = file_path.with_suffix(".txt.tmp")

        # When
        atomic_write(file_path, "content")

        # Then
        assert file_path.exists()
        assert not temp_path.exists()

    def test_handles_unicode_content(self, tmp_path: Path) -> None:
        """[P2] Given unicode content, when atomic_write called, then content preserved."""
        # Given
        file_path = tmp_path / "test.txt"
        content = "中文内容 with emoji 🎉"

        # When
        atomic_write(file_path, content)

        # Then
        assert file_path.read_text(encoding="utf-8") == content


# =============================================================================
# Test atomic_write_json
# =============================================================================

class TestAtomicWriteJson:
    """Tests for atomic_write_json function."""

    def test_writes_json_to_file(self, tmp_path: Path) -> None:
        """[P1] Given dict data, when atomic_write_json called, then valid JSON written."""
        # Given
        file_path = tmp_path / "test.json"
        data = {"key": "value", "number": 42}

        # When
        atomic_write_json(file_path, data)

        # Then
        assert file_path.exists()
        loaded = json.loads(file_path.read_text())
        assert loaded == data

    def test_preserves_unicode_in_json(self, tmp_path: Path) -> None:
        """[P2] Given unicode in dict, when written, then unicode preserved."""
        # Given
        file_path = tmp_path / "test.json"
        data = {"message": "中文消息", "emoji": "🚀"}

        # When
        atomic_write_json(file_path, data)

        # Then
        loaded = json.loads(file_path.read_text(encoding="utf-8"))
        assert loaded == data

    def test_uses_specified_indent(self, tmp_path: Path) -> None:
        """[P2] Given indent parameter, when written, then JSON is indented."""
        # Given
        file_path = tmp_path / "test.json"
        data = {"key": "value"}

        # When
        atomic_write_json(file_path, data, indent=4)

        # Then
        content = file_path.read_text()
        assert "    " in content  # 4-space indent present

    def test_handles_nested_data(self, tmp_path: Path) -> None:
        """[P2] Given nested dict, when written, then structure preserved."""
        # Given
        file_path = tmp_path / "test.json"
        data = {
            "outer": {
                "inner": {
                    "deep": ["a", "b", "c"]
                }
            }
        }

        # When
        atomic_write_json(file_path, data)

        # Then
        loaded = json.loads(file_path.read_text())
        assert loaded == data

    def test_file_ends_with_newline(self, tmp_path: Path) -> None:
        """[P2] Given any data, when written, then file ends with newline."""
        # Given
        file_path = tmp_path / "test.json"
        data = {"key": "value"}

        # When
        atomic_write_json(file_path, data)

        # Then
        content = file_path.read_text()
        assert content.endswith("\n")


# =============================================================================
# Test atomic_append
# =============================================================================

class TestAtomicAppend:
    """Tests for atomic_append function."""

    def test_creates_file_if_not_exists(self, tmp_path: Path) -> None:
        """[P1] Given non-existent file, when atomic_append called, then file created."""
        # Given
        file_path = tmp_path / "test.txt"
        content = "First line"

        # When
        atomic_append(file_path, content)

        # Then
        assert file_path.exists()
        assert file_path.read_text() == content

    def test_appends_to_existing_file(self, tmp_path: Path) -> None:
        """[P1] Given existing file, when atomic_append called, then content appended."""
        # Given
        file_path = tmp_path / "test.txt"
        file_path.write_text("First\n")
        new_content = "Second\n"

        # When
        atomic_append(file_path, new_content)

        # Then
        assert file_path.read_text() == "First\nSecond\n"

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        """[P2] Given nested path, when atomic_append called, then directories created."""
        # Given
        file_path = tmp_path / "subdir" / "test.txt"

        # When
        atomic_append(file_path, "content")

        # Then
        assert file_path.exists()

    def test_handles_multiple_appends(self, tmp_path: Path) -> None:
        """[P2] Given multiple appends, when called sequentially, then all content present."""
        # Given
        file_path = tmp_path / "test.txt"
        lines = ["Line 1\n", "Line 2\n", "Line 3\n"]

        # When
        for line in lines:
            atomic_append(file_path, line)

        # Then
        assert file_path.read_text() == "".join(lines)
