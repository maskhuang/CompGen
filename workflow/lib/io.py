"""
CompGene I/O Utilities.

This module provides atomic file writing and other I/O utilities
for the CompGene pipeline.

Key features:
- Atomic writes: Ensures checkpoint safety for Snakemake --rerun-incomplete
- Temporary file cleanup: Automatically cleans up .tmp files on failure

Source: ADR-004 Log Format (Atomic Write Pattern)
"""

import json
import logging
from pathlib import Path
from typing import Any, Optional

# Module logger
_logger = logging.getLogger(__name__)


def atomic_write(path: Path, content: str, encoding: str = "utf-8") -> None:
    """
    Atomically write text content to a file.

    Writes to a temporary file first, then renames to the target path.
    This ensures that readers never see a partially written file, which is
    critical for Snakemake's --rerun-incomplete to work correctly.

    On failure, the temporary file is automatically cleaned up to prevent
    leaving orphaned .tmp files.

    Args:
        path: Target file path.
        content: Text content to write.
        encoding: Text encoding (default: utf-8).

    Raises:
        OSError: If write or rename fails (after cleanup).

    Example:
        >>> atomic_write(Path("output.txt"), "Hello, World!")
    """
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        temp_path.write_text(content, encoding=encoding)
        temp_path.rename(path)
    except Exception:
        # Clean up temporary file on failure
        temp_path.unlink(missing_ok=True)
        raise


def atomic_write_json(path: Path, data: Any, indent: int = 2) -> None:
    """
    Atomically write JSON data to a file.

    Args:
        path: Target file path.
        data: JSON-serializable data.
        indent: JSON indentation level (default: 2).

    Example:
        >>> atomic_write_json(Path("data.json"), {"key": "value"})
    """
    content = json.dumps(data, ensure_ascii=False, indent=indent)
    atomic_write(path, content + "\n")


def atomic_write_bytes(path: Path, content: bytes) -> None:
    """
    Atomically write binary content to a file.

    Writes to a temporary file first, then renames to the target path.
    This ensures that readers never see a partially written file.

    On failure, the temporary file is automatically cleaned up.

    Args:
        path: Target file path.
        content: Binary content to write.

    Raises:
        OSError: If write or rename fails (after cleanup).

    Example:
        >>> atomic_write_bytes(Path("output.bin"), b"\\x00\\x01\\x02")
    """
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        temp_path.write_bytes(content)
        temp_path.rename(path)
    except Exception:
        # Clean up temporary file on failure
        temp_path.unlink(missing_ok=True)
        raise


def atomic_append(path: Path, content: str, encoding: str = "utf-8") -> None:
    """
    Append content to a file, creating it if it doesn't exist.

    Note: This is NOT atomic for the append operation itself,
    but ensures the file is created atomically if new.

    Args:
        path: Target file path.
        content: Text content to append.
        encoding: Text encoding (default: utf-8).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding=encoding) as f:
        f.write(content)


def cleanup_temp_files(
    directory: Path,
    pattern: str = "*.tmp",
    recursive: bool = True
) -> list[Path]:
    """
    Clean up orphaned temporary files in a directory.

    This function removes .tmp files that may have been left behind
    after interrupted operations, ensuring clean state for
    Snakemake --rerun-incomplete.

    Args:
        directory: Directory to clean.
        pattern: Glob pattern for temp files (default: "*.tmp").
        recursive: Whether to search recursively (default: True).

    Returns:
        List of paths that were removed.

    Example:
        >>> removed = cleanup_temp_files(Path("results/"))
        >>> print(f"Cleaned up {len(removed)} temp files")
    """
    if not directory.exists():
        return []

    glob_method = directory.rglob if recursive else directory.glob
    removed: list[Path] = []

    for temp_file in glob_method(pattern):
        if temp_file.is_file():
            try:
                temp_file.unlink()
                removed.append(temp_file)
            except OSError as e:
                # Log but don't fail - cleanup is best-effort
                _logger.debug("Could not remove temp file %s: %s", temp_file, e)

    return removed


def get_temp_path(path: Path) -> Path:
    """
    Get the temporary file path for a given target path.

    This follows the convention used by atomic_write functions.

    Args:
        path: Target file path.

    Returns:
        Path to the corresponding .tmp file.

    Example:
        >>> get_temp_path(Path("results/output.tsv"))
        Path('results/output.tsv.tmp')
    """
    return path.with_suffix(path.suffix + ".tmp")
