"""
CompGene Checksum Utilities.

This module provides checksum computation for files and directories.
Used for audit tracking and cache invalidation.

Source: ADR-004 Audit Granularity
"""

import hashlib
from pathlib import Path
from typing import Optional

# Default chunk size for large file processing (8 MB)
DEFAULT_CHUNK_SIZE = 8 * 1024 * 1024


def compute_file_checksum(
    path: Path,
    algorithm: str = "sha256",
    chunk_size: int = DEFAULT_CHUNK_SIZE
) -> str:
    """
    Compute checksum for a single file.

    Uses chunked reading to handle large files without memory issues.

    Args:
        path: Path to the file.
        algorithm: Hash algorithm name (default: sha256).
        chunk_size: Size of chunks to read (default: 8 MB).

    Returns:
        Checksum string in format "{algorithm}:{hex_digest}".

    Raises:
        FileNotFoundError: If the file doesn't exist.
        IsADirectoryError: If the path is a directory.

    Example:
        >>> checksum = compute_file_checksum(Path("genome.fa"))
        >>> print(checksum)
        sha256:abc123...
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    if path.is_dir():
        raise IsADirectoryError(f"Path is a directory, use compute_dir_checksum: {path}")

    hasher = hashlib.new(algorithm)

    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            hasher.update(chunk)

    return f"{algorithm}:{hasher.hexdigest()}"


def compute_dir_checksum(
    path: Path,
    algorithm: str = "sha256",
    pattern: str = "*",
    recursive: bool = True
) -> str:
    """
    Compute a combined checksum for all files in a directory.

    Files are sorted by name to ensure deterministic results.
    The checksum is computed over the concatenation of:
    - Relative file path (UTF-8 encoded)
    - File content

    Args:
        path: Path to the directory.
        algorithm: Hash algorithm name (default: sha256).
        pattern: Glob pattern to match files (default: "*").
        recursive: Whether to search recursively (default: True).

    Returns:
        Checksum string in format "{algorithm}:{hex_digest}".

    Raises:
        FileNotFoundError: If the directory doesn't exist.
        NotADirectoryError: If the path is not a directory.

    Example:
        >>> checksum = compute_dir_checksum(Path("proteins/"))
        >>> print(checksum)
        sha256:def456...
    """
    if not path.exists():
        raise FileNotFoundError(f"Directory not found: {path}")
    if not path.is_dir():
        raise NotADirectoryError(f"Path is not a directory: {path}")

    hasher = hashlib.new(algorithm)

    # Get files matching pattern, sorted for deterministic order
    glob_method = path.rglob if recursive else path.glob
    files = sorted(f for f in glob_method(pattern) if f.is_file())

    for file_path in files:
        # Include relative path in hash for structure tracking
        rel_path = file_path.relative_to(path)
        hasher.update(str(rel_path).encode("utf-8"))

        # Include file content
        with file_path.open("rb") as f:
            while True:
                chunk = f.read(DEFAULT_CHUNK_SIZE)
                if not chunk:
                    break
                hasher.update(chunk)

    return f"{algorithm}:{hasher.hexdigest()}"


def compute_checksums(
    paths: dict[str, Path],
    algorithm: str = "sha256"
) -> dict[str, str]:
    """
    Compute checksums for multiple files/directories.

    Automatically detects whether each path is a file or directory
    and uses the appropriate checksum function.

    Args:
        paths: Dictionary mapping names to paths.
        algorithm: Hash algorithm name (default: sha256).

    Returns:
        Dictionary mapping names to checksum strings.

    Example:
        >>> checksums = compute_checksums({
        ...     "genome": Path("genome.fa"),
        ...     "proteins": Path("proteins/")
        ... })
    """
    result: dict[str, str] = {}

    for name, path in paths.items():
        if not path.exists():
            result[name] = f"{algorithm}:MISSING"
        elif path.is_dir():
            result[name] = compute_dir_checksum(path, algorithm)
        else:
            result[name] = compute_file_checksum(path, algorithm)

    return result
