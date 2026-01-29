"""
CompGene Checksum Utilities.

This module provides checksum computation for files and directories.
Used for audit tracking and cache invalidation.

Features:
- Efficient chunked reading for large files
- Optional checksum caching for fast re-runs (< 1 sec cache hit)
- Directory checksum with deterministic ordering

Source: ADR-004 Audit Granularity
"""

import hashlib
import json
import logging
from pathlib import Path
from typing import Optional

# Module logger
_logger = logging.getLogger(__name__)

# Default chunk size for large file processing (8 MB)
DEFAULT_CHUNK_SIZE = 8 * 1024 * 1024

# Default cache file name
DEFAULT_CACHE_FILE = ".checksum_cache.json"


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


# =============================================================================
# Checksum Cache (Story 1.7 - Performance Optimization)
# =============================================================================

class ChecksumCache:
    """
    A file-based cache for checksums to speed up re-runs.

    Stores checksums keyed by (path, mtime) to avoid recomputing
    checksums for unchanged files. This enables < 1 second cache
    hit performance as required by NFR4.

    The cache is stored as a JSON file in the specified cache directory.

    Thread Safety:
        This class is NOT thread-safe. In Snakemake parallel execution,
        each rule should use its own ChecksumCache instance or ensure
        exclusive access through external synchronization. The save()
        method uses atomic writes to prevent corruption from concurrent
        saves, but the in-memory cache is not protected.

    Attributes:
        cache_file: Path to the cache JSON file.
        cache: In-memory cache dictionary.
        dirty: Whether the cache has unsaved changes.

    Example:
        >>> cache = ChecksumCache(Path("results/meta"))
        >>> checksum = cache.get_or_compute(Path("genome.fa"))
        >>> cache.save()  # Persist to disk
    """

    def __init__(
        self,
        cache_dir: Path,
        cache_file_name: str = DEFAULT_CACHE_FILE
    ):
        """
        Initialize the checksum cache.

        Args:
            cache_dir: Directory to store the cache file.
            cache_file_name: Name of the cache file.
        """
        self.cache_file = cache_dir / cache_file_name
        self.cache: dict[str, dict[str, str | float]] = {}
        self.dirty = False
        self._load()

    def _load(self) -> None:
        """Load cache from disk if it exists."""
        if self.cache_file.exists():
            try:
                with self.cache_file.open("r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except (json.JSONDecodeError, OSError):
                # Corrupted or unreadable cache - start fresh
                self.cache = {}

    def save(self) -> None:
        """Save cache to disk if there are unsaved changes."""
        if not self.dirty:
            return

        self.cache_file.parent.mkdir(parents=True, exist_ok=True)
        temp_path = self.cache_file.with_suffix(".json.tmp")
        try:
            with temp_path.open("w", encoding="utf-8") as f:
                json.dump(self.cache, f, ensure_ascii=False, indent=2)
            temp_path.rename(self.cache_file)
            self.dirty = False
        except OSError:
            temp_path.unlink(missing_ok=True)
            raise

    def _get_cache_key(self, path: Path) -> str:
        """Generate a cache key from the absolute path."""
        return str(path.resolve())

    def _get_mtime(self, path: Path) -> float:
        """Get the modification time of a file."""
        return path.stat().st_mtime

    def get(self, path: Path) -> Optional[str]:
        """
        Get a cached checksum if still valid.

        Args:
            path: Path to the file.

        Returns:
            Cached checksum if valid, None otherwise.
        """
        if not path.exists() or path.is_dir():
            return None

        key = self._get_cache_key(path)
        entry = self.cache.get(key)

        if entry is None:
            return None

        # Check if mtime matches (file unchanged)
        try:
            current_mtime = self._get_mtime(path)
            cached_mtime = entry.get("mtime", 0)

            if current_mtime == cached_mtime:
                checksum = entry.get("checksum")
                if isinstance(checksum, str):
                    return checksum
                return None
        except OSError as e:
            _logger.debug("Could not get mtime for %s: %s", path, e)

        return None

    def set(self, path: Path, checksum: str) -> None:
        """
        Store a checksum in the cache.

        Args:
            path: Path to the file.
            checksum: The computed checksum.
        """
        if not path.exists() or path.is_dir():
            return

        key = self._get_cache_key(path)
        try:
            mtime = self._get_mtime(path)
            self.cache[key] = {
                "checksum": checksum,
                "mtime": mtime
            }
            self.dirty = True
        except OSError as e:
            _logger.debug("Could not cache checksum for %s: %s", path, e)

    def get_or_compute(
        self,
        path: Path,
        algorithm: str = "sha256",
        chunk_size: int = DEFAULT_CHUNK_SIZE
    ) -> str:
        """
        Get checksum from cache or compute if not cached/stale.

        This is the primary method for efficient checksum operations.

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
            >>> cache = ChecksumCache(Path("results/meta"))
            >>> cs1 = cache.get_or_compute(Path("genome.fa"))  # Computes
            >>> cs2 = cache.get_or_compute(Path("genome.fa"))  # Cache hit
        """
        # Try cache first
        cached = self.get(path)
        if cached is not None:
            return cached

        # Compute checksum
        checksum = compute_file_checksum(path, algorithm, chunk_size)

        # Store in cache
        self.set(path, checksum)

        return checksum

    def invalidate(self, path: Path) -> bool:
        """
        Remove a path from the cache.

        Args:
            path: Path to invalidate.

        Returns:
            True if entry was removed, False if not in cache.
        """
        key = self._get_cache_key(path)
        if key in self.cache:
            del self.cache[key]
            self.dirty = True
            return True
        return False

    def clear(self) -> None:
        """Clear all cached entries."""
        if self.cache:
            self.cache = {}
            self.dirty = True

    def __len__(self) -> int:
        """Return number of cached entries."""
        return len(self.cache)


def compute_checksums_cached(
    paths: dict[str, Path],
    cache: ChecksumCache,
    algorithm: str = "sha256"
) -> dict[str, str]:
    """
    Compute checksums for multiple files using cache.

    Uses the provided cache for efficient re-computation.
    Directories are not cached (computed each time).

    Args:
        paths: Dictionary mapping names to paths.
        cache: ChecksumCache instance to use.
        algorithm: Hash algorithm name (default: sha256).

    Returns:
        Dictionary mapping names to checksum strings.

    Example:
        >>> cache = ChecksumCache(Path("results/meta"))
        >>> checksums = compute_checksums_cached({
        ...     "genome": Path("genome.fa"),
        ...     "proteins": Path("proteins/")
        ... }, cache)
        >>> cache.save()
    """
    result: dict[str, str] = {}

    for name, path in paths.items():
        if not path.exists():
            result[name] = f"{algorithm}:MISSING"
        elif path.is_dir():
            # Directories are not cached
            result[name] = compute_dir_checksum(path, algorithm)
        else:
            result[name] = cache.get_or_compute(path, algorithm)

    return result
