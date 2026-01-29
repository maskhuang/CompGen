"""
CompGene I/O Utilities.

This module provides atomic file writing and other I/O utilities
for the CompGene pipeline.

Source: ADR-004 Log Format (Atomic Write Pattern)
"""

import json
from pathlib import Path
from typing import Any


def atomic_write(path: Path, content: str, encoding: str = "utf-8") -> None:
    """
    Atomically write text content to a file.

    Writes to a temporary file first, then renames to the target path.
    This ensures that readers never see a partially written file.

    Args:
        path: Target file path.
        content: Text content to write.
        encoding: Text encoding (default: utf-8).

    Example:
        >>> atomic_write(Path("output.txt"), "Hello, World!")
    """
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path.write_text(content, encoding=encoding)
    temp_path.rename(path)


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
