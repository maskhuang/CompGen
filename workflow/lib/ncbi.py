"""
CompGene NCBI Download Module.

This module provides functionality for downloading genomic data from NCBI,
including genome sequences, annotations, and protein sequences.

Features:
- NCBI accession parsing (GCF_/GCA_ formats)
- Rate-limited API access (3 req/s default, configurable with API key)
- Exponential backoff retry logic
- Local caching to avoid redundant downloads
- Atomic file writes for checkpoint safety

Source: Story 2.3 - NCBI Annotation Download
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import threading
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Optional, TypeVar

from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.io import atomic_write_json

# Module logger
_logger = logging.getLogger(__name__)

# Type variable for retry function
T = TypeVar("T")

# =============================================================================
# Constants
# =============================================================================

# NCBI FTP base URL for genomes
NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# NCBI Datasets API base URL
NCBI_DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# NCBI Entrez E-utilities base URL
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Default rate limit (requests per second) without API key
DEFAULT_RATE_LIMIT = 3.0

# Rate limit with API key (10 req/s)
API_KEY_RATE_LIMIT = 10.0

# Default cache manifest filename
CACHE_MANIFEST = "download_manifest.json"

# User agent for NCBI requests
USER_AGENT = "CompGene/1.0 (Comparative Genomics Pipeline)"

# =============================================================================
# Accession Types and Parsing
# =============================================================================


class AccessionType(str, Enum):
    """NCBI accession type enum."""

    REFSEQ = "refseq"    # GCF_ prefix - RefSeq assemblies
    GENBANK = "genbank"  # GCA_ prefix - GenBank assemblies


@dataclass
class NCBIAccession:
    """
    Parsed NCBI assembly accession.

    Represents a parsed NCBI assembly accession with components for
    constructing FTP paths and API queries.

    Attributes:
        raw: Original input string (may include NCBI: prefix)
        prefix: GCF or GCA
        numeric: 9-digit numeric identifier
        version: Optional version number
        accession_type: RefSeq or GenBank

    Example:
        >>> acc = parse_ncbi_accession("NCBI:GCF_000165445.1")
        >>> acc.full_accession
        'GCF_000165445.1'
        >>> acc.ftp_path_parts
        ('000', '165', '445')
    """

    raw: str
    prefix: str
    numeric: str
    version: Optional[str]
    accession_type: AccessionType

    @property
    def full_accession(self) -> str:
        """Return the full accession without NCBI: prefix."""
        if self.version:
            return f"{self.prefix}_{self.numeric}.{self.version}"
        return f"{self.prefix}_{self.numeric}"

    @property
    def base_accession(self) -> str:
        """Return the accession without version."""
        return f"{self.prefix}_{self.numeric}"

    @property
    def ftp_path_parts(self) -> tuple[str, str, str]:
        """
        Return the three-part path segments for NCBI FTP.

        NCBI FTP organizes files by splitting the numeric ID into
        three 3-digit segments.

        Returns:
            Tuple of (part1, part2, part3) for FTP path construction.

        Example:
            >>> acc = parse_ncbi_accession("GCF_000165445")
            >>> acc.ftp_path_parts
            ('000', '165', '445')
        """
        return (
            self.numeric[0:3],
            self.numeric[3:6],
            self.numeric[6:9]
        )

    @property
    def ftp_base_path(self) -> str:
        """
        Return the base FTP path for this accession.

        Returns:
            FTP path up to the accession directory level.

        Example:
            >>> acc = parse_ncbi_accession("GCF_000165445.1")
            >>> acc.ftp_base_path
            'GCF/000/165/445'
        """
        p1, p2, p3 = self.ftp_path_parts
        return f"{self.prefix}/{p1}/{p2}/{p3}"


# Regex pattern for NCBI accession parsing
# Matches: NCBI:GCF_000165445, GCF_000165445.1, GCA_000165445, etc.
NCBI_PATTERN = re.compile(
    r"^(?:NCBI:)?(GC[FA])_(\d{9})(?:\.(\d+))?$",
    re.IGNORECASE
)


def parse_ncbi_accession(identifier: str) -> NCBIAccession:
    """
    Parse an NCBI accession identifier.

    Accepts accessions in the following formats:
    - NCBI:GCF_000165445.1 (with prefix)
    - GCF_000165445.1 (RefSeq with version)
    - GCA_000165445 (GenBank without version)

    Args:
        identifier: The accession string to parse.

    Returns:
        Parsed NCBIAccession object.

    Raises:
        CompGeneError: If the format is invalid (E_INPUT_FORMAT).

    Example:
        >>> acc = parse_ncbi_accession("NCBI:GCF_000165445.1")
        >>> acc.prefix
        'GCF'
        >>> acc.version
        '1'
    """
    if not identifier:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            "Empty NCBI accession",
            details="Accession identifier cannot be empty"
        )

    match = NCBI_PATTERN.match(identifier.strip())
    if not match:
        raise CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            f"Invalid NCBI accession format: {identifier}",
            details=(
                "Expected format: NCBI:GCF_XXXXXXXXX.V or GCA_XXXXXXXXX "
                "where X is a digit and V is the version"
            )
        )

    prefix = match.group(1).upper()
    numeric = match.group(2)
    version = match.group(3)

    acc_type = AccessionType.REFSEQ if prefix == "GCF" else AccessionType.GENBANK

    return NCBIAccession(
        raw=identifier,
        prefix=prefix,
        numeric=numeric,
        version=version,
        accession_type=acc_type
    )


def is_ncbi_accession(identifier: Optional[str]) -> bool:
    """
    Check if a string is a valid NCBI accession.

    Args:
        identifier: String to check (can be None).

    Returns:
        True if the string matches NCBI accession format.
    """
    if not identifier:
        return False
    return NCBI_PATTERN.match(identifier.strip()) is not None


# =============================================================================
# Rate Limiter
# =============================================================================


class RateLimiter:
    """
    Token bucket rate limiter for NCBI API compliance.

    Implements a thread-safe token bucket algorithm to enforce
    rate limits on API requests. Default rate is 3 requests/second
    as per NCBI guidelines, or 10 req/s with an API key.

    Attributes:
        rate: Maximum requests per second.
        tokens: Current available tokens.
        lock: Threading lock for thread safety.

    Example:
        >>> limiter = RateLimiter(3.0)
        >>> limiter.acquire()  # Returns immediately if tokens available
        >>> limiter.acquire()  # May wait if rate limit reached
    """

    def __init__(self, requests_per_second: float = DEFAULT_RATE_LIMIT):
        """
        Initialize the rate limiter.

        Args:
            requests_per_second: Maximum requests per second (default: 3.0).
        """
        self.rate = requests_per_second
        self.tokens = requests_per_second
        self.last_update = time.monotonic()
        self.lock = threading.Lock()

    def acquire(self) -> float:
        """
        Acquire a token, blocking if necessary.

        Blocks until a token is available, ensuring the rate limit
        is not exceeded.

        Returns:
            The time spent waiting (0.0 if no wait was needed).
        """
        with self.lock:
            now = time.monotonic()
            elapsed = now - self.last_update

            # Refill tokens based on elapsed time
            self.tokens = min(self.rate, self.tokens + elapsed * self.rate)
            self.last_update = now

            if self.tokens >= 1:
                self.tokens -= 1
                return 0.0

            # Calculate wait time needed
            wait_time = (1 - self.tokens) / self.rate
            time.sleep(wait_time)

            self.tokens = 0
            self.last_update = time.monotonic()
            return wait_time

    def set_rate(self, requests_per_second: float) -> None:
        """
        Update the rate limit.

        Args:
            requests_per_second: New rate limit.
        """
        with self.lock:
            self.rate = requests_per_second


# Global rate limiter instance
_rate_limiter: Optional[RateLimiter] = None
_rate_limiter_lock = threading.Lock()


def get_rate_limiter(
    api_key: Optional[str] = None,
    rate_limit: Optional[float] = None
) -> RateLimiter:
    """
    Get or create the global rate limiter.

    Args:
        api_key: Optional NCBI API key to increase rate limit.
        rate_limit: Optional custom rate limit (overrides api_key-based detection).

    Returns:
        RateLimiter instance.
    """
    global _rate_limiter

    with _rate_limiter_lock:
        if _rate_limiter is None:
            if rate_limit is not None:
                rate = rate_limit
            else:
                rate = API_KEY_RATE_LIMIT if api_key else DEFAULT_RATE_LIMIT
            _rate_limiter = RateLimiter(rate)
        elif rate_limit is not None:
            # Update rate if custom rate provided
            _rate_limiter.set_rate(rate_limit)
        elif api_key:
            # Upgrade rate if API key provided
            _rate_limiter.set_rate(API_KEY_RATE_LIMIT)

        return _rate_limiter


def reset_rate_limiter() -> None:
    """
    Reset the global rate limiter.

    Used primarily for testing to ensure clean state between tests.
    """
    global _rate_limiter

    with _rate_limiter_lock:
        _rate_limiter = None


# =============================================================================
# Retry Logic
# =============================================================================


class RetryableError(Exception):
    """Exception wrapper for retryable errors."""
    pass


def with_retry(
    func: Callable[[], T],
    max_retries: int = 3,
    base_delay: float = 1.0,
    max_delay: float = 60.0,
    retryable_exceptions: tuple = (
        urllib.error.URLError,
        TimeoutError,
        ConnectionError,
        RetryableError,
    )
) -> T:
    """
    Execute a function with exponential backoff retry.

    Retries the function on specified exceptions using exponential
    backoff strategy: delay = min(base_delay * 2^attempt, max_delay).

    Args:
        func: Function to execute (no arguments).
        max_retries: Maximum number of retry attempts (default: 3).
        base_delay: Initial delay in seconds (default: 1.0).
        max_delay: Maximum delay in seconds (default: 60.0).
        retryable_exceptions: Tuple of exception types to retry.

    Returns:
        The function's return value.

    Raises:
        The last exception if all retries fail.

    Example:
        >>> def fetch_data():
        ...     return urllib.request.urlopen(url).read()
        >>> data = with_retry(fetch_data, max_retries=3)
    """
    last_exception: Optional[Exception] = None

    for attempt in range(max_retries + 1):
        try:
            return func()
        except retryable_exceptions as e:
            last_exception = e

            if attempt < max_retries:
                delay = min(base_delay * (2 ** attempt), max_delay)
                _logger.warning(
                    f"Retry {attempt + 1}/{max_retries} after {delay:.1f}s: {e}"
                )
                time.sleep(delay)

    # All retries exhausted
    if last_exception is not None:
        raise last_exception
    raise RuntimeError("Unexpected retry state")


# =============================================================================
# NCBI API Key Support
# =============================================================================


def get_ncbi_api_key() -> Optional[str]:
    """
    Get NCBI API key from environment.

    Checks the NCBI_API_KEY environment variable.

    Returns:
        API key string or None if not set.
    """
    return os.environ.get("NCBI_API_KEY")


def build_url_with_api_key(url: str, api_key: Optional[str] = None) -> str:
    """
    Add API key to URL query parameters if available.

    Args:
        url: Base URL.
        api_key: Optional API key (uses env var if None).

    Returns:
        URL with api_key parameter if key is available.
    """
    key = api_key or get_ncbi_api_key()
    if not key:
        return url

    separator = "&" if "?" in url else "?"
    return f"{url}{separator}api_key={key}"


# =============================================================================
# Download Functions
# =============================================================================


@dataclass
class DownloadResult:
    """
    Result of a file download operation.

    Attributes:
        url: Source URL.
        path: Local file path.
        size_bytes: File size in bytes.
        sha256: SHA256 checksum of the file.
        content_type: HTTP content type (if available).
    """

    url: str
    path: str
    size_bytes: int
    sha256: str
    content_type: Optional[str] = None


def download_file(
    url: str,
    dest_path: Path,
    rate_limiter: Optional[RateLimiter] = None,
    chunk_size: int = 8192,
    timeout: int = 300
) -> DownloadResult:
    """
    Download a file from URL to local path.

    Uses atomic write pattern (temp file + rename) to ensure
    checkpoint safety. Computes SHA256 checksum during download.

    Args:
        url: Source URL.
        dest_path: Destination file path.
        rate_limiter: Optional rate limiter instance.
        chunk_size: Download chunk size (default: 8192).
        timeout: Request timeout in seconds (default: 300).

    Returns:
        DownloadResult with file metadata.

    Raises:
        CompGeneError: On download failure.
    """
    if rate_limiter:
        rate_limiter.acquire()

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = dest_path.with_suffix(dest_path.suffix + ".tmp")

    sha256 = hashlib.sha256()
    total_bytes = 0
    content_type: Optional[str] = None

    try:
        request = urllib.request.Request(url)
        request.add_header("User-Agent", USER_AGENT)

        with urllib.request.urlopen(request, timeout=timeout) as response:
            content_type = response.headers.get("Content-Type")

            with open(temp_path, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    sha256.update(chunk)
                    total_bytes += len(chunk)

        # Atomic rename
        temp_path.rename(dest_path)

        _logger.debug(f"Downloaded {total_bytes} bytes to {dest_path}")

        return DownloadResult(
            url=url,
            path=str(dest_path),
            size_bytes=total_bytes,
            sha256=sha256.hexdigest(),
            content_type=content_type
        )

    except urllib.error.HTTPError as e:
        temp_path.unlink(missing_ok=True)
        if e.code == 429:  # Rate limit
            raise CompGeneError(
                ErrorCode.E_NET_RATE_LIMIT,
                f"NCBI rate limit exceeded: {url}",
                details=str(e)
            )
        elif e.code == 404:
            raise CompGeneError(
                ErrorCode.E_INPUT_MISSING,
                f"File not found on NCBI: {url}",
                details=str(e)
            )
        else:
            raise CompGeneError(
                ErrorCode.E_NONZERO_EXIT,
                f"HTTP error downloading {url}: {e.code}",
                details=str(e)
            )

    except urllib.error.URLError as e:
        temp_path.unlink(missing_ok=True)
        raise CompGeneError(
            ErrorCode.E_TIMEOUT,
            f"Network error downloading {url}",
            details=str(e)
        )

    except Exception as e:
        temp_path.unlink(missing_ok=True)
        raise


def fetch_url_content(
    url: str,
    rate_limiter: Optional[RateLimiter] = None,
    timeout: int = 60
) -> bytes:
    """
    Fetch URL content as bytes.

    Args:
        url: URL to fetch.
        rate_limiter: Optional rate limiter.
        timeout: Request timeout in seconds.

    Returns:
        Response content as bytes.

    Raises:
        urllib.error.URLError: On network error.
        urllib.error.HTTPError: On HTTP error.
    """
    if rate_limiter:
        rate_limiter.acquire()

    request = urllib.request.Request(url)
    request.add_header("User-Agent", USER_AGENT)

    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


# =============================================================================
# NCBI FTP Directory Resolution
# =============================================================================


def resolve_ftp_directory(
    accession: NCBIAccession,
    rate_limiter: Optional[RateLimiter] = None
) -> str:
    """
    Resolve the full FTP directory path for an accession.

    NCBI FTP directories include the assembly name after the accession,
    so we need to list the directory to find the exact path.

    Args:
        accession: Parsed NCBI accession.
        rate_limiter: Optional rate limiter.

    Returns:
        Full FTP directory URL.

    Raises:
        CompGeneError: If the directory cannot be resolved.
    """
    base_url = f"{NCBI_FTP_BASE}/{accession.ftp_base_path}"

    def _fetch():
        content = fetch_url_content(base_url + "/", rate_limiter, timeout=60)
        return content.decode("utf-8")

    try:
        html_content = with_retry(_fetch, max_retries=3)
    except Exception as e:
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Cannot access NCBI FTP directory for {accession.full_accession}",
            details=str(e)
        )

    # Parse HTML directory listing
    # Look for links matching the accession pattern
    pattern = re.compile(
        rf'href="({re.escape(accession.full_accession)}[^"]*)"',
        re.IGNORECASE
    )

    matches = pattern.findall(html_content)

    if not matches:
        # Try base accession without version
        pattern_base = re.compile(
            rf'href="({re.escape(accession.base_accession)}[^"]*)"',
            re.IGNORECASE
        )
        matches = pattern_base.findall(html_content)

    if not matches:
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"No assembly found for {accession.full_accession}",
            details=f"Directory listing at {base_url} does not contain matching assembly"
        )

    # Use the first match (should be the directory)
    dir_name = matches[0].rstrip("/")
    return f"{base_url}/{dir_name}"


def get_assembly_files(
    ftp_dir: str,
    rate_limiter: Optional[RateLimiter] = None
) -> dict[str, str]:
    """
    Get available files in an NCBI assembly directory.

    Args:
        ftp_dir: FTP directory URL.
        rate_limiter: Optional rate limiter.

    Returns:
        Dictionary mapping file type to URL.
        Keys: 'genome', 'annotation', 'protein', 'report'
    """
    def _fetch():
        content = fetch_url_content(ftp_dir + "/", rate_limiter, timeout=60)
        return content.decode("utf-8")

    try:
        html_content = with_retry(_fetch, max_retries=3)
    except Exception as e:
        raise CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            f"Cannot list files in {ftp_dir}",
            details=str(e)
        )

    # Extract file links
    link_pattern = re.compile(r'href="([^"]+)"')
    files = link_pattern.findall(html_content)

    result: dict[str, str] = {}

    for filename in files:
        filename = filename.rstrip("/")

        if filename.endswith("_genomic.fna.gz"):
            # Exclude masked genome files and RNA/CDS extracts
            if "_rm." not in filename and "_sm." not in filename and \
               "_rna_from_" not in filename and "_cds_from_" not in filename:
                result["genome"] = f"{ftp_dir}/{filename}"

        elif filename.endswith("_genomic.gff.gz"):
            result["annotation"] = f"{ftp_dir}/{filename}"

        elif filename.endswith("_protein.faa.gz"):
            result["protein"] = f"{ftp_dir}/{filename}"

        elif filename.endswith("_assembly_report.txt"):
            result["report"] = f"{ftp_dir}/{filename}"

    return result


# =============================================================================
# Cache Management
# =============================================================================


def is_cached(accession: str, cache_dir: Path) -> bool:
    """
    Check if an accession's data is cached.

    Checks for the existence of the download manifest and all
    referenced files.

    Args:
        accession: Full accession string (without NCBI: prefix).
        cache_dir: Base cache directory.

    Returns:
        True if all cached files exist.
    """
    manifest_path = cache_dir / accession / CACHE_MANIFEST

    if not manifest_path.exists():
        return False

    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        # Check all files exist
        for file_info in manifest.get("files", []):
            file_path = Path(file_info.get("path", ""))
            if not file_path.exists():
                _logger.debug(f"Cached file missing: {file_path}")
                return False

        return True

    except (json.JSONDecodeError, KeyError, OSError) as e:
        _logger.debug(f"Cache check failed: {e}")
        return False


def verify_cache(accession: str, cache_dir: Path) -> bool:
    """
    Verify cached files using checksums.

    Args:
        accession: Full accession string.
        cache_dir: Base cache directory.

    Returns:
        True if all checksums match.
    """
    manifest_path = cache_dir / accession / CACHE_MANIFEST

    if not manifest_path.exists():
        return False

    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        for file_info in manifest.get("files", []):
            file_path = Path(file_info.get("path", ""))
            expected_sha256 = file_info.get("sha256")

            if not file_path.exists():
                return False

            if expected_sha256:
                # Compute actual checksum
                sha256 = hashlib.sha256()
                with open(file_path, "rb") as f:
                    while True:
                        chunk = f.read(8192)
                        if not chunk:
                            break
                        sha256.update(chunk)

                if sha256.hexdigest() != expected_sha256:
                    _logger.warning(
                        f"Checksum mismatch for {file_path}: "
                        f"expected {expected_sha256}, got {sha256.hexdigest()}"
                    )
                    return False

        return True

    except Exception as e:
        _logger.debug(f"Cache verification failed: {e}")
        return False


def write_download_manifest(
    accession: str,
    output_dir: Path,
    files: list[dict[str, Any]],
    metadata: Optional[dict[str, Any]] = None
) -> Path:
    """
    Write a download manifest file.

    Args:
        accession: Full accession string.
        output_dir: Directory to write manifest.
        files: List of file info dictionaries.
        metadata: Optional additional metadata.

    Returns:
        Path to the manifest file.
    """
    manifest = {
        "accession": accession,
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "compgene_version": "1.0.0",
        "files": files,
        "metadata": metadata or {}
    }

    manifest_path = output_dir / CACHE_MANIFEST
    atomic_write_json(manifest_path, manifest)

    return manifest_path


# =============================================================================
# Main Download Function
# =============================================================================


@dataclass
class NCBIDownloadResult:
    """
    Result of downloading an NCBI dataset.

    Attributes:
        accession: Full accession string.
        output_dir: Directory containing downloaded files.
        genome_path: Path to genome file (if downloaded).
        annotation_path: Path to annotation file (if downloaded).
        protein_path: Path to protein file (if downloaded).
        manifest_path: Path to download manifest.
        from_cache: Whether files were retrieved from cache.
    """

    accession: str
    output_dir: Path
    genome_path: Optional[Path] = None
    annotation_path: Optional[Path] = None
    protein_path: Optional[Path] = None
    manifest_path: Optional[Path] = None
    from_cache: bool = False


def download_ncbi_dataset(
    accession: str,
    output_dir: Path,
    force: bool = False,
    api_key: Optional[str] = None,
    rate_limit: Optional[float] = None,
    download_genome: bool = True,
    download_annotation: bool = True,
    download_protein: bool = True
) -> NCBIDownloadResult:
    """
    Download a complete NCBI dataset for an accession.

    Downloads genome, annotation, and optionally protein files from
    NCBI FTP. Uses caching to avoid redundant downloads.

    Args:
        accession: NCBI accession (with or without NCBI: prefix).
        output_dir: Base output directory.
        force: Force re-download even if cached (default: False).
        api_key: Optional NCBI API key for higher rate limits.
        rate_limit: Optional custom rate limit in requests per second.
        download_genome: Download genome file (default: True).
        download_annotation: Download annotation file (default: True).
        download_protein: Download protein file (default: True).

    Returns:
        NCBIDownloadResult with paths to downloaded files.

    Raises:
        CompGeneError: On parsing, network, or download errors.

    Example:
        >>> result = download_ncbi_dataset(
        ...     "NCBI:GCF_000165445.1",
        ...     Path("downloads"),
        ...     force=False
        ... )
        >>> print(result.genome_path)
    """
    # Parse accession
    parsed = parse_ncbi_accession(accession)
    full_accession = parsed.full_accession

    # Set up output directory
    accession_dir = output_dir / full_accession
    accession_dir.mkdir(parents=True, exist_ok=True)

    # Check cache
    if not force and is_cached(full_accession, output_dir):
        _logger.info(f"Using cached data for {full_accession}")

        # Load manifest to get paths
        manifest_path = accession_dir / CACHE_MANIFEST
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)

        result = NCBIDownloadResult(
            accession=full_accession,
            output_dir=accession_dir,
            manifest_path=manifest_path,
            from_cache=True
        )

        for file_info in manifest.get("files", []):
            file_type = file_info.get("type")
            file_path = Path(file_info.get("path", ""))

            if file_type == "genome":
                result.genome_path = file_path
            elif file_type == "annotation":
                result.annotation_path = file_path
            elif file_type == "protein":
                result.protein_path = file_path

        return result

    # Get rate limiter
    rate_limiter = get_rate_limiter(api_key, rate_limit)

    _logger.info(f"Resolving FTP directory for {full_accession}")

    # Resolve FTP directory
    ftp_dir = resolve_ftp_directory(parsed, rate_limiter)
    _logger.info(f"FTP directory: {ftp_dir}")

    # Get available files
    available_files = get_assembly_files(ftp_dir, rate_limiter)
    _logger.info(f"Available files: {list(available_files.keys())}")

    # Download files
    downloaded_files: list[dict[str, Any]] = []
    result = NCBIDownloadResult(
        accession=full_accession,
        output_dir=accession_dir
    )

    # Download genome
    if download_genome and "genome" in available_files:
        genome_url = available_files["genome"]
        genome_dest = accession_dir / "genome.fna.gz"

        _logger.info(f"Downloading genome: {genome_url}")

        def _download_genome():
            return download_file(genome_url, genome_dest, rate_limiter)

        dl_result = with_retry(_download_genome, max_retries=3)

        downloaded_files.append({
            "type": "genome",
            "url": dl_result.url,
            "path": str(genome_dest),
            "size_bytes": dl_result.size_bytes,
            "sha256": dl_result.sha256
        })
        result.genome_path = genome_dest

    # Download annotation
    if download_annotation and "annotation" in available_files:
        annotation_url = available_files["annotation"]
        annotation_dest = accession_dir / "annotation.gff.gz"

        _logger.info(f"Downloading annotation: {annotation_url}")

        def _download_annotation():
            return download_file(annotation_url, annotation_dest, rate_limiter)

        dl_result = with_retry(_download_annotation, max_retries=3)

        downloaded_files.append({
            "type": "annotation",
            "url": dl_result.url,
            "path": str(annotation_dest),
            "size_bytes": dl_result.size_bytes,
            "sha256": dl_result.sha256
        })
        result.annotation_path = annotation_dest

    # Download protein
    if download_protein and "protein" in available_files:
        protein_url = available_files["protein"]
        protein_dest = accession_dir / "protein.faa.gz"

        _logger.info(f"Downloading protein: {protein_url}")

        def _download_protein():
            return download_file(protein_url, protein_dest, rate_limiter)

        dl_result = with_retry(_download_protein, max_retries=3)

        downloaded_files.append({
            "type": "protein",
            "url": dl_result.url,
            "path": str(protein_dest),
            "size_bytes": dl_result.size_bytes,
            "sha256": dl_result.sha256
        })
        result.protein_path = protein_dest

    # Write manifest
    manifest_path = write_download_manifest(
        full_accession,
        accession_dir,
        downloaded_files,
        metadata={
            "ftp_directory": ftp_dir,
            "accession_type": parsed.accession_type.value
        }
    )
    result.manifest_path = manifest_path

    _logger.info(f"Download complete: {len(downloaded_files)} files")

    return result


# =============================================================================
# NCBIClient Class (Higher-level Interface)
# =============================================================================


class NCBIClient:
    """
    High-level client for NCBI data access.

    Provides a convenient interface for downloading and managing
    NCBI genomic data with caching and rate limiting.

    Attributes:
        cache_dir: Directory for caching downloads.
        api_key: Optional NCBI API key.
        rate_limiter: Rate limiter instance.

    Example:
        >>> client = NCBIClient(Path("downloads"))
        >>> result = client.download("NCBI:GCF_000165445.1")
        >>> print(result.annotation_path)
    """

    def __init__(
        self,
        cache_dir: Path,
        api_key: Optional[str] = None,
        rate_limit: Optional[float] = None
    ):
        """
        Initialize NCBI client.

        Args:
            cache_dir: Directory for caching downloads.
            api_key: Optional NCBI API key (default: from env).
            rate_limit: Custom rate limit (default: auto-detect).
        """
        self.cache_dir = Path(cache_dir)
        self.api_key = api_key or get_ncbi_api_key()

        if rate_limit:
            self.rate_limiter = RateLimiter(rate_limit)
        else:
            rate = API_KEY_RATE_LIMIT if self.api_key else DEFAULT_RATE_LIMIT
            self.rate_limiter = RateLimiter(rate)

    def download(
        self,
        accession: str,
        force: bool = False,
        download_genome: bool = True,
        download_annotation: bool = True,
        download_protein: bool = True
    ) -> NCBIDownloadResult:
        """
        Download an NCBI dataset.

        Args:
            accession: NCBI accession string.
            force: Force re-download (default: False).
            download_genome: Download genome file (default: True).
            download_annotation: Download annotation file (default: True).
            download_protein: Download protein file (default: True).

        Returns:
            NCBIDownloadResult with file paths.
        """
        return download_ncbi_dataset(
            accession=accession,
            output_dir=self.cache_dir,
            force=force,
            api_key=self.api_key,
            download_genome=download_genome,
            download_annotation=download_annotation,
            download_protein=download_protein
        )

    def is_cached(self, accession: str) -> bool:
        """
        Check if an accession is cached.

        Args:
            accession: NCBI accession string.

        Returns:
            True if cached.
        """
        parsed = parse_ncbi_accession(accession)
        return is_cached(parsed.full_accession, self.cache_dir)

    def verify_cache(self, accession: str) -> bool:
        """
        Verify cached data integrity.

        Args:
            accession: NCBI accession string.

        Returns:
            True if cache is valid.
        """
        parsed = parse_ncbi_accession(accession)
        return verify_cache(parsed.full_accession, self.cache_dir)

    def clear_cache(self, accession: str) -> bool:
        """
        Clear cached data for an accession.

        Args:
            accession: NCBI accession string.

        Returns:
            True if cache was cleared.
        """
        parsed = parse_ncbi_accession(accession)
        cache_path = self.cache_dir / parsed.full_accession

        if cache_path.exists():
            shutil.rmtree(cache_path)
            return True
        return False
