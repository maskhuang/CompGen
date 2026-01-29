"""
Unit tests for workflow/lib/checksum.py

Tests checksum computation for files and directories.
"""

import tempfile
from pathlib import Path

import pytest

from workflow.lib.checksum import (
    compute_file_checksum,
    compute_dir_checksum,
    compute_checksums,
)


# =============================================================================
# Test compute_file_checksum
# =============================================================================

class TestComputeFileChecksum:
    """Tests for compute_file_checksum function."""

    def test_computes_sha256_by_default(self, tmp_path: Path) -> None:
        """[P1] Given a file, when compute_file_checksum called, then sha256 returned."""
        # Given
        file_path = tmp_path / "test.txt"
        file_path.write_text("Hello, World!")

        # When
        checksum = compute_file_checksum(file_path)

        # Then
        assert checksum.startswith("sha256:")
        assert len(checksum.split(":")[1]) == 64  # SHA256 hex length

    def test_returns_consistent_checksum(self, tmp_path: Path) -> None:
        """[P1] Given same content, when computed twice, then same checksum returned."""
        # Given
        file_path = tmp_path / "test.txt"
        file_path.write_text("Consistent content")

        # When
        checksum1 = compute_file_checksum(file_path)
        checksum2 = compute_file_checksum(file_path)

        # Then
        assert checksum1 == checksum2

    def test_different_content_different_checksum(self, tmp_path: Path) -> None:
        """[P1] Given different content, when computed, then different checksums."""
        # Given
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_text("Content A")
        file2.write_text("Content B")

        # When
        checksum1 = compute_file_checksum(file1)
        checksum2 = compute_file_checksum(file2)

        # Then
        assert checksum1 != checksum2

    def test_raises_for_nonexistent_file(self, tmp_path: Path) -> None:
        """[P1] Given nonexistent file, when called, then FileNotFoundError raised."""
        # Given
        file_path = tmp_path / "nonexistent.txt"

        # When/Then
        with pytest.raises(FileNotFoundError):
            compute_file_checksum(file_path)

    def test_raises_for_directory(self, tmp_path: Path) -> None:
        """[P1] Given directory path, when called, then IsADirectoryError raised."""
        # Given
        dir_path = tmp_path / "subdir"
        dir_path.mkdir()

        # When/Then
        with pytest.raises(IsADirectoryError):
            compute_file_checksum(dir_path)

    def test_supports_md5_algorithm(self, tmp_path: Path) -> None:
        """[P2] Given md5 algorithm, when called, then md5 checksum returned."""
        # Given
        file_path = tmp_path / "test.txt"
        file_path.write_text("Test content")

        # When
        checksum = compute_file_checksum(file_path, algorithm="md5")

        # Then
        assert checksum.startswith("md5:")
        assert len(checksum.split(":")[1]) == 32  # MD5 hex length

    def test_handles_large_file_in_chunks(self, tmp_path: Path) -> None:
        """[P2] Given large file, when computed, then checksum returned without memory issue."""
        # Given - Create a file larger than default chunk size
        file_path = tmp_path / "large.bin"
        # Write 1 MB of data
        file_path.write_bytes(b"x" * (1024 * 1024))

        # When
        checksum = compute_file_checksum(file_path, chunk_size=1024)  # Small chunk

        # Then
        assert checksum.startswith("sha256:")

    def test_handles_empty_file(self, tmp_path: Path) -> None:
        """[P2] Given empty file, when computed, then valid checksum returned."""
        # Given
        file_path = tmp_path / "empty.txt"
        file_path.write_text("")

        # When
        checksum = compute_file_checksum(file_path)

        # Then
        assert checksum.startswith("sha256:")
        # SHA256 of empty string
        assert "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" in checksum


# =============================================================================
# Test compute_dir_checksum
# =============================================================================

class TestComputeDirChecksum:
    """Tests for compute_dir_checksum function."""

    def test_computes_checksum_for_directory(self, tmp_path: Path) -> None:
        """[P1] Given directory with files, when called, then sha256 returned."""
        # Given
        dir_path = tmp_path / "testdir"
        dir_path.mkdir()
        (dir_path / "file1.txt").write_text("Content 1")
        (dir_path / "file2.txt").write_text("Content 2")

        # When
        checksum = compute_dir_checksum(dir_path)

        # Then
        assert checksum.startswith("sha256:")

    def test_consistent_for_same_content(self, tmp_path: Path) -> None:
        """[P1] Given same directory content, when computed twice, then same checksum."""
        # Given
        dir_path = tmp_path / "testdir"
        dir_path.mkdir()
        (dir_path / "file.txt").write_text("Content")

        # When
        checksum1 = compute_dir_checksum(dir_path)
        checksum2 = compute_dir_checksum(dir_path)

        # Then
        assert checksum1 == checksum2

    def test_different_for_different_content(self, tmp_path: Path) -> None:
        """[P1] Given different directory content, then different checksums."""
        # Given
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()
        (dir1 / "file.txt").write_text("Content A")
        (dir2 / "file.txt").write_text("Content B")

        # When
        checksum1 = compute_dir_checksum(dir1)
        checksum2 = compute_dir_checksum(dir2)

        # Then
        assert checksum1 != checksum2

    def test_raises_for_nonexistent_directory(self, tmp_path: Path) -> None:
        """[P1] Given nonexistent directory, when called, then FileNotFoundError."""
        # Given
        dir_path = tmp_path / "nonexistent"

        # When/Then
        with pytest.raises(FileNotFoundError):
            compute_dir_checksum(dir_path)

    def test_raises_for_file_path(self, tmp_path: Path) -> None:
        """[P1] Given file path, when called, then NotADirectoryError."""
        # Given
        file_path = tmp_path / "file.txt"
        file_path.write_text("content")

        # When/Then
        with pytest.raises(NotADirectoryError):
            compute_dir_checksum(file_path)

    def test_includes_subdirectories_by_default(self, tmp_path: Path) -> None:
        """[P2] Given nested directories, when recursive=True, then all files included."""
        # Given
        dir_path = tmp_path / "testdir"
        subdir = dir_path / "subdir"
        subdir.mkdir(parents=True)
        (dir_path / "file1.txt").write_text("Root")
        (subdir / "file2.txt").write_text("Nested")

        # When
        checksum_recursive = compute_dir_checksum(dir_path, recursive=True)
        checksum_flat = compute_dir_checksum(dir_path, recursive=False)

        # Then
        assert checksum_recursive != checksum_flat

    def test_handles_empty_directory(self, tmp_path: Path) -> None:
        """[P2] Given empty directory, when computed, then valid checksum returned."""
        # Given
        dir_path = tmp_path / "empty"
        dir_path.mkdir()

        # When
        checksum = compute_dir_checksum(dir_path)

        # Then
        assert checksum.startswith("sha256:")

    def test_file_order_is_deterministic(self, tmp_path: Path) -> None:
        """[P2] Given files created in different order, then same checksum."""
        # Given - Create two directories with same files in different order
        dir1 = tmp_path / "dir1"
        dir2 = tmp_path / "dir2"
        dir1.mkdir()
        dir2.mkdir()

        # Create files in different order
        (dir1 / "aaa.txt").write_text("A")
        (dir1 / "zzz.txt").write_text("Z")

        (dir2 / "zzz.txt").write_text("Z")
        (dir2 / "aaa.txt").write_text("A")

        # When
        checksum1 = compute_dir_checksum(dir1)
        checksum2 = compute_dir_checksum(dir2)

        # Then - Should be same because files are sorted
        assert checksum1 == checksum2


# =============================================================================
# Test compute_checksums
# =============================================================================

class TestComputeChecksums:
    """Tests for compute_checksums function."""

    def test_computes_multiple_checksums(self, tmp_path: Path) -> None:
        """[P1] Given multiple paths, when called, then all checksums returned."""
        # Given
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.write_text("Content 1")
        file2.write_text("Content 2")

        paths = {"first": file1, "second": file2}

        # When
        checksums = compute_checksums(paths)

        # Then
        assert "first" in checksums
        assert "second" in checksums
        assert checksums["first"].startswith("sha256:")
        assert checksums["second"].startswith("sha256:")

    def test_handles_mixed_files_and_directories(self, tmp_path: Path) -> None:
        """[P1] Given mix of files and dirs, when called, then appropriate method used."""
        # Given
        file_path = tmp_path / "file.txt"
        dir_path = tmp_path / "dir"
        file_path.write_text("File content")
        dir_path.mkdir()
        (dir_path / "nested.txt").write_text("Nested")

        paths = {"file": file_path, "directory": dir_path}

        # When
        checksums = compute_checksums(paths)

        # Then
        assert checksums["file"].startswith("sha256:")
        assert checksums["directory"].startswith("sha256:")

    def test_marks_missing_as_missing(self, tmp_path: Path) -> None:
        """[P1] Given nonexistent path, when called, then marked as MISSING."""
        # Given
        missing_path = tmp_path / "nonexistent.txt"
        existing_path = tmp_path / "exists.txt"
        existing_path.write_text("content")

        paths = {"missing": missing_path, "exists": existing_path}

        # When
        checksums = compute_checksums(paths)

        # Then
        assert "MISSING" in checksums["missing"]
        assert "MISSING" not in checksums["exists"]

    def test_returns_empty_dict_for_empty_input(self) -> None:
        """[P2] Given empty paths dict, when called, then empty dict returned."""
        # Given
        paths: dict[str, Path] = {}

        # When
        checksums = compute_checksums(paths)

        # Then
        assert checksums == {}
